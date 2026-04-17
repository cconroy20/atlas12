! ==============================================================================
! mod_mklinelist.f90 — Kurucz line list preprocessor, embedded module
!
! Replaces the standalone mklinelist program and the legacy setup.sh pipeline:
!   synbeg | rgfall | rpredict | rmolecasc (×N) | rh2ofast
!
! Public interface:
!   call run_mklinelist(wlbeg, wlend, resolu, lines_list_path)
!
! After the call, the module-level arrays lte_lines(:) and nlte_lines(:) are
! allocated and populated.  nlines_lte and nlines_nlte hold their sizes.
! These arrays replace the fort.12 and fort.19 binary scratch files.
!
! Line list sources are specified in a plain-text file (lines.list by default):
!   gfall   /path/to/gfall.dat
!   predict /path/to/predict.bin
!   mol     /path/to/molecule1.dat
!   mol     /path/to/molecule2.dat
!   h2o     /path/to/h2o.bin
! Blank lines and lines beginning with # are ignored.
! Order of entries is irrelevant; readers always run in canonical order:
!   gfall -> predict -> mol(s) -> h2o
! The file must exist; a hard stop is issued if it is not found.
! ==============================================================================

module mod_mklinelist
  implicit none
  private

  ! ============================================================================
  !  PUBLIC TYPES
  ! ============================================================================

  !  LTE line record — mirrors the fort.12 binary record layout:
  !    NBUFF(I4)  CGF(R4)  NELION(I4)  ELO(R4)
  !    GAMRF(R4)  GAMSF(R4)  GAMWF(R4)          [28 bytes]
  type, public :: lte_line_t
    integer  :: nbuff    ! grid-index offset from ixwlbeg
    real(4)  :: cgf      ! line-strength factor (units: per (cm^-1 g^-1 cm^2))
    integer  :: nelion   ! species index (nelem-1)*6 + ion_stage
    real(4)  :: elo      ! lower level energy (cm^-1)
    real(4)  :: gamrf    ! radiative  damping / 4πν
    real(4)  :: gamsf    ! Stark      damping / 4πν
    real(4)  :: gamwf    ! van der Waals damping / 4πν
  end type lte_line_t

  !  NLTE / complex-profile line record — mirrors the fort.19 binary record:
  !    WLVAC(R8)  ELO(R4)  GF(R4)  NBLO(I4)  NBUP(I4)
  !    NELION(I4)  TYPE(I4)  NCON(I4)  NELIONX(I4)
  !    GAMMAR(R4)  GAMMAS(R4)  GAMMAW(R4)  NBUFF(I4)  LIM(I4)  [60 bytes]
  type, public :: nlte_line_t
    real(8)  :: wlvac    ! vacuum wavelength (nm)
    real(4)  :: elo      ! lower level energy (cm^-1)
    real(4)  :: gf       ! oscillator strength (CGF for most types; raw GF for
                         !   autoionizing TYPE=1 and continuum TYPE>3)
    integer  :: nblo     ! lower NLTE level index
    integer  :: nbup     ! upper NLTE level index
    integer  :: nelion   ! species index
    integer  :: itype    ! line type flag (see read_gfall dispatch table)
    integer  :: ncon     ! number of continuum edges (isotope field reuse)
    integer  :: nelionx  ! species index for departure coefficients
    real(4)  :: gammar   ! radiative  damping / 4πν
    real(4)  :: gammas   ! Stark      damping / 4πν
    real(4)  :: gammaw   ! van der Waals damping / 4πν
    integer  :: nbuff    ! grid-index offset from ixwlbeg
    integer  :: lim      ! wing extent index (0=widest, 6=narrowest)
  end type nlte_line_t

  ! ============================================================================
  !  PUBLIC MODULE ARRAYS AND COUNTERS
  ! ============================================================================

  type(lte_line_t),  allocatable, public :: lte_lines(:)
  type(nlte_line_t), allocatable, public :: nlte_lines(:)
  integer,                        public :: nlines_lte  = 0
  integer,                        public :: nlines_nlte = 0

  ! ============================================================================
  !  PUBLIC ENTRY POINT
  ! ============================================================================

  public :: run_mklinelist

  ! ============================================================================
  !  MODULE-LEVEL IONISATION POTENTIAL TABLE  (replaces COMMON /potion/)
  ! ============================================================================

  real(8), save :: potion(999)

  ! ============================================================================
  !  TEMPERATURE GATE FOR COOL-ATMOSPHERE MOLECULAR LINES
  !
  !  Line lists for species that only form in cool atmospheres (H2O, TiO) are
  !  skipped when Teff exceeds this threshold, even if they are listed in
  !  lines.list.  This avoids wasted I/O and memory on molecular opacities
  !  that contribute negligibly in warmer stars.
  ! ============================================================================

  real(8), parameter :: TEFF_COOL_LIMIT = 5000.0d0

contains

  ! ============================================================================
  !  RUN_MKLINELIST — top-level entry point called from SYNTHE
  !
  !  Reads lines.list, dispatches to readers in canonical order, assembles
  !  the lte_lines and nlte_lines module arrays.
  ! ============================================================================
  subroutine run_mklinelist(wlbeg, wlend, resolu, teff, lines_list_path, datadir)
    real(8),          intent(in) :: wlbeg, wlend, resolu, teff
    character(len=*), intent(in) :: lines_list_path
    character(len=*), intent(in) :: datadir

    ! synbeg_params — grid scalars passed to every reader
    real(8) :: ratio, ratiolg
    integer :: length, ixwlbeg
    real(8) :: wbegin

    ! File paths parsed from lines.list
    character(len=512) :: gfall_file, predict_file, h2o_file
    character(len=512), save :: mol_files(256)
    integer :: nmol

    ! Temporary per-reader arrays (LTE)
    type(lte_line_t),  allocatable :: lte_gfall(:),   lte_predict(:)
    type(lte_line_t),  allocatable :: lte_mol(:),     lte_h2o(:)
    integer :: n_lte_gfall, n_lte_predict, n_lte_mol, n_lte_h2o

    ! Temporary per-reader arrays (NLTE)
    type(nlte_line_t), allocatable :: nlte_gfall(:)
    integer :: n_nlte_gfall

    integer :: imol, ioff
    integer :: ixwlend
    integer :: n_lte_mol_before   ! per-file delta for mol entries

    ! ---- initialise ionisation potential table ----
    call ionpots()

    ! ---- compute derived grid quantities ----
    ratio   = 1.0d0 + 1.0d0 / resolu
    ratiolg = log(ratio)

    ixwlbeg = int(log(wlbeg) / ratiolg)
    wbegin  = exp(ixwlbeg * ratiolg)
    if (wbegin < wlbeg) then
      ixwlbeg = ixwlbeg + 1
      wbegin  = exp(ixwlbeg * ratiolg)
    end if

    ixwlend = int(log(wlend) / ratiolg)
    if (exp(ixwlend * ratiolg) >= wlend) ixwlend = ixwlend - 1
    length  = ixwlend - ixwlbeg + 1

    ! ---- parse lines.list ----
    gfall_file   = ''
    predict_file = ''
    h2o_file     = ''
    nmol         = 0
    call parse_lines_list(lines_list_path, datadir, &
                          gfall_file, predict_file, h2o_file, &
                          mol_files, nmol)

    ! ---- initialise per-reader output arrays ----
    n_lte_gfall   = 0;  allocate(lte_gfall(0))
    n_lte_predict = 0;  allocate(lte_predict(0))
    n_lte_mol     = 0;  allocate(lte_mol(0))
    n_lte_h2o     = 0;  allocate(lte_h2o(0))
    n_nlte_gfall  = 0;  allocate(nlte_gfall(0))

    ! ---- header ----
    write(6,'(/a)') '=== mklinelist ==='
    write(6,'(a4,2x,a12,2x,a12,2x,a)') 'src', 'n_lte', 'n_nlte', 'file'

    ! ---- canonical order: gfall -> predict -> mol(s) -> h2o ----

    if (gfall_file /= '') then
      call read_gfall(gfall_file, wlbeg, wlend, ratiolg, ixwlbeg, &
                      lte_gfall, n_lte_gfall, nlte_gfall, n_nlte_gfall)
      write(6,'(a4,2x,i12,2x,i12,2x,a)') 'gf', n_lte_gfall, n_nlte_gfall, &
        trim(basename(gfall_file))
    end if

    if (predict_file /= '') then
      call read_predict(predict_file, wlbeg, wlend, ratiolg, ixwlbeg, &
                        lte_predict, n_lte_predict)
      write(6,'(a4,2x,i12,2x,i12,2x,a)') 'pr', n_lte_predict, 0, &
        trim(basename(predict_file))
    end if

    do imol = 1, nmol
      ! Skip TiO files when Teff is above the cool-atmosphere limit.
      if (is_tio_file(mol_files(imol)) .and. teff > TEFF_COOL_LIMIT) then
        write(6,'(a4,2x,a12,2x,a12,2x,a,a)') 'mol', 'skip', '', &
          'Teff > 5000 K: ', trim(basename(mol_files(imol)))
        cycle
      end if
      n_lte_mol_before = n_lte_mol
      call read_molec(mol_files(imol), wlbeg, wlend, ratiolg, ixwlbeg, &
                      lte_mol, n_lte_mol)
      write(6,'(a4,2x,i12,2x,i12,2x,a)') 'mol', n_lte_mol - n_lte_mol_before, 0, &
        trim(basename(mol_files(imol)))
    end do

    if (h2o_file /= '' .and. teff > TEFF_COOL_LIMIT) then
      write(6,'(a4,2x,a12,2x,a12,2x,a,a)') 'h2o', 'skip', '', &
        'Teff > 5000 K: ', trim(basename(h2o_file))
    else if (h2o_file /= '') then
      call read_h2o(h2o_file, wlbeg, wlend, ratiolg, ixwlbeg, &
                    lte_h2o, n_lte_h2o)
      write(6,'(a4,2x,i12,2x,i12,2x,a)') 'h2o', n_lte_h2o, 0, &
        trim(basename(h2o_file))
    end if

    ! ---- assemble module arrays from per-reader temporaries ----
    nlines_nlte = n_nlte_gfall
    nlines_lte  = n_lte_gfall + n_lte_predict + n_lte_mol + n_lte_h2o

    allocate(nlte_lines(nlines_nlte))
    if (nlines_nlte > 0) &
      nlte_lines(1:nlines_nlte) = nlte_gfall(1:nlines_nlte)

    allocate(lte_lines(nlines_lte))
    ioff = 0
    if (n_lte_gfall > 0) then
      lte_lines(ioff+1 : ioff+n_lte_gfall) = lte_gfall(1:n_lte_gfall)
      ioff = ioff + n_lte_gfall
    end if
    if (n_lte_predict > 0) then
      lte_lines(ioff+1 : ioff+n_lte_predict) = lte_predict(1:n_lte_predict)
      ioff = ioff + n_lte_predict
    end if
    if (n_lte_mol > 0) then
      lte_lines(ioff+1 : ioff+n_lte_mol) = lte_mol(1:n_lte_mol)
      ioff = ioff + n_lte_mol
    end if
    if (n_lte_h2o > 0) then
      lte_lines(ioff+1 : ioff+n_lte_h2o) = lte_h2o(1:n_lte_h2o)
    end if

    deallocate(lte_gfall, lte_predict, lte_mol, lte_h2o, nlte_gfall)

    write(6,'(a4,2x,a12,2x,a12)') '---', '------------', '------------'
    write(6,'(a4,2x,i12,2x,i12)') 'tot', nlines_lte, nlines_nlte
    write(6,'(a)') ''

  end subroutine run_mklinelist


  ! ============================================================================
  !  PARSE_LINES_LIST — read lines.list and populate file path variables
  ! ============================================================================
  subroutine parse_lines_list(listfile, datadir, gfall_file, predict_file, h2o_file, &
                               mol_files, nmol)
    character(len=*),   intent(in)  :: listfile
    character(len=*),   intent(in)  :: datadir
    character(len=512), intent(out) :: gfall_file, predict_file, h2o_file
    character(len=512), intent(out) :: mol_files(256)
    integer,            intent(out) :: nmol

    integer, parameter :: LU = 50
    character(len=1024) :: line
    character(len=32)   :: keyword
    character(len=512)  :: filepath
    integer :: ios

    gfall_file   = ''
    predict_file = ''
    h2o_file     = ''
    nmol         = 0

    open(unit=LU, file=listfile, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') 'ERROR: lines.list not found: ', trim(listfile)
      stop 1
    end if

    do
      read(LU, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (line == '')            cycle
      if (line(1:1) == '#')      cycle

      ! Read keyword only via list-directed (stops at first whitespace).
      ! List-directed read of the full line would truncate filepath at '/'
      ! since Fortran treats '/' as a record terminator in list-directed I/O.
      read(line, *, iostat=ios) keyword
      if (ios /= 0) then
        write(6,'(a,a)') 'WARNING: ignoring malformed line in lines.list: ', trim(line)
        cycle
      end if
      ! Extract filepath as everything after the keyword, stripping any
      ! leading whitespace (spaces or tabs) between keyword and filename.
      ! adjustl only strips spaces, so we scan manually.
      block
        integer :: i, klen
        klen = len_trim(keyword)
        i    = klen + 1
        do while (i <= len_trim(line))
          if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
            i = i + 1
          else
            exit
          end if
        end do
        filepath = trim(datadir) // line(i:len_trim(line))
      end block

      select case (trim(keyword))
      case ('gfall')
        gfall_file = trim(filepath)
      case ('predict')
        predict_file = trim(filepath)
      case ('mol')
        if (nmol < 256) then
          nmol = nmol + 1
          mol_files(nmol) = trim(filepath)
        else
          write(6,'(a)') 'WARNING: more than 256 mol entries in lines.list; extras ignored'
        end if
      case ('h2o')
        h2o_file = trim(filepath)
      case default
        write(6,'(a,a)') 'WARNING: unrecognised keyword in lines.list: ', trim(keyword)
      end select
    end do

    close(LU)

  end subroutine parse_lines_list


  ! ============================================================================
  !  READ_GFALL — read Kurucz gfall ASCII atomic/ionic line list
  !               (translated from rgfall.for)
  !
  !  LTE plain Voigt lines go into lte_arr.
  !  Complex-profile lines (H, He, autoionizing, PRD, continuum edges, NLTE
  !  departure coefficient lines) go into nlte_arr.
  ! ============================================================================
  subroutine read_gfall(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                         lte_arr, n_lte, nlte_arr, n_nlte)
    character(len=*),              intent(in)    :: filename
    real(8),                       intent(in)    :: wlbeg, wlend, ratiolg
    integer,                       intent(in)    :: ixwlbeg
    type(lte_line_t),  allocatable, intent(inout) :: lte_arr(:)
    integer,                       intent(inout) :: n_lte
    type(nlte_line_t), allocatable, intent(inout) :: nlte_arr(:)
    integer,                       intent(inout) :: n_nlte

    ! ---- LINDAT fields ----
    real(8)  :: wl, e, ep, label(2), labelp(2), other1(2), other2(2)
    real(8)  :: wlvac
    real(4)  :: nelion_r
    real(4)  :: ref, gflog, xj, xjp, code, gf, gs, gr, gw
    real(4)  :: dwl, dgflog, dgammar, dgammas, dgammaw, dwliso
    integer  :: nblo, nbup, iso1, iso2, isoshift
    real(4)  :: x1, x2

    real(8)  :: lindat8(14)
    real(4)  :: lindat4(28)
    equivalence (lindat8(1), wl),      (lindat4(1), nelion_r)

    character(len=10) :: cother1, cother2
    equivalence (cother1, other1(1)), (cother2, other2(1))
    character(len=3)  :: auto_flag
    character(len=6)  :: ixfixfp

    ! CODEX: species eligible for NLTE departure coefficients
    real(4), parameter :: CODEX(17) = [ &
      1.0, 2.0, 2.01, 6.0, 6.01, 12.0, 12.01, 13.0, 13.01, &
      14.0, 14.01, 20.0, 20.01, 8.0, 11.0, 5.0, 19.0]

    real(8), parameter :: DELLIM(7) = &
      [100.d0, 30.d0, 10.d0, 3.d0, 1.d0, 0.3d0, 0.1d0]

    ! ---- working variables ----
    integer  :: ios, nelem, icharge, linesize, lim, ncon, nelionx, itype, ic
    integer  :: ixwl, nbuff, ishift, ishiftp
    real(8)  :: delfactor, eshift, eshiftp, frelin, cgf, frq4pi
    real(8)  :: effnsq, zeff, rsqup, rsqlo, eup
    real(8)  :: gammar_d, gammas_d, gammaw_d, elo_d
    integer  :: nelion_i

    ! ---- dynamic growth buffers ----
    integer, parameter :: CHUNK = 100000
    type(lte_line_t),  allocatable :: lte_buf(:)
    type(nlte_line_t), allocatable :: nlte_buf(:)
    integer :: lte_cap, nlte_cap

    lte_cap  = CHUNK;  allocate(lte_buf(lte_cap))
    nlte_cap = CHUNK;  allocate(nlte_buf(nlte_cap))
    n_lte  = 0
    n_nlte = 0

    open(unit=11, file=filename, status='old', action='read', &
         form='formatted', iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') 'ERROR: cannot open gfall file: ', trim(filename); stop 1
    end if

    delfactor = 1.0d0
    if (wlbeg > 500.0d0) delfactor = wlbeg / 500.0d0

    dwl = 0.0;  dgflog = 0.0;  dgammar = 0.0
    dgammas = 0.0;  dgammaw = 0.0;  dwliso = 0.0
    other1(2) = 0.0d0;  other2(1) = 0.0d0;  other2(2) = 0.0d0

    do
      read(11, &
        '(F11.4,F7.3,F6.2,F12.3,F5.1,1X,A8,A2,F12.3,F5.1,1X,A8,A2,' // &
        'F6.2,F6.2,F6.2,A4,I2,I2,I3,F6.3,I3,F6.3,A8,A2,A8,A2,2I5,I6)', &
        iostat=ios) &
        wl, gflog, code, e, xj, label, ep, xjp, labelp, &
        gr, gs, gw, ref, nblo, nbup, iso1, x1, iso2, x2, &
        other1, other2, lim, lim, isoshift
      if (ios < 0) exit
      if (ios > 0) then
        write(6,'(a)') ' WARNING: read_gfall: read error, stopping line read'; exit
      end if

      read(cother1, '(2I5)', iostat=ios) ishift, ishiftp
      if (ios /= 0) then; ishift = 0;  ishiftp = 0; end if
      read(cother2, '(A6,I1,A3)', iostat=ios) ixfixfp, linesize, auto_flag
      if (ios /= 0) then; linesize = 0; auto_flag = '   '; end if

      eshift  = ishift  * 0.001d0
      eshiftp = ishiftp * 0.001d0
      dwliso  = real(isoshift, 4) * 0.0001     ! milli-Angstrom -> nm

      wlvac = abs(wl) + dwl
      if (transfer(labelp(1), ' ') == 'CONTINUU') &
        wlvac = 1.d7 / abs(abs(ep) + eshiftp - (abs(e) + eshift)) + dwl + dwliso

      ! Guard: skip any record with non-physical wlvac (e.g. continuum edge
      ! where the energy level difference is zero giving wlvac=Inf or NaN)
      if (.not. (wlvac > 0.0d0 .and. wlvac < 1.0d8)) cycle

      if (wlvac > wlend + dellim(1)) exit

      ixwl  = int(log(wlvac) / ratiolg + 0.5d0)
      nbuff = ixwl - ixwlbeg + 1

      lim = min(8 - linesize, 7)
      if (code == 1.0) lim = 1

      if (wlvac < wlbeg - dellim(lim) * delfactor) cycle
      if (wlvac > wlend + dellim(lim) * delfactor) cycle
      if (auto_flag == 'COR') cycle

      ! ---- oscillator strength and energy ----
      gf    = 10.0**(gflog + dgflog + x1 + x2)
      elo_d = dble(min(abs(e), abs(ep)))

      ! ---- damping constants ----
      gammar_d = 10.0d0**(gr + dgammar)
      gammas_d = 10.0d0**(gs + dgammas)
      gammaw_d = 10.0d0**(gw + dgammaw)
      if (auto_flag == 'AUT' .and. gs > 0.0) &
        gammas_d = -10.0d0**(-gs + dgammas)

      if (gr == 0.0) gammar_d = 2.223d13 / wlvac**2

      ! ---- NELION ----
      nelem    = int(code)
      icharge  = nint((code - real(nelem,4)) * 100.0)
      zeff     = icharge + 1
      nelion_i = nelem*6 - 6 + int(zeff)
      if (nelem > 19 .and. nelem < 29 .and. icharge > 5) &
        nelion_i = 6*(nelem + icharge*10 - 30) - 1

      ! ---- default Stark width ----
      if (gs == 0.0) then
        if (code < 100.0) then
          eup    = dble(max(abs(e), abs(ep)))
          effnsq = 25.0d0
          call ionpot_index(nelem, icharge, eup, effnsq, zeff)
          gammas_d = 1.0d-8 * effnsq**2 * sqrt(effnsq)
          gs = real(log10(gammas_d), 4)
        else
          gammas_d = 1.0d-5
        end if
      end if

      ! ---- default van der Waals width ----
      if (gw == 0.0) then
        if (code < 100.0) then
          eup    = dble(max(abs(e), abs(ep)))
          effnsq = 25.0d0
          call ionpot_index(nelem, icharge, eup, effnsq, zeff)
          effnsq = min(effnsq, 1000.0d0)
          rsqup  = 2.5d0 * (effnsq / zeff)**2
          call ionpot_index_lo(nelem, icharge, elo_d, effnsq, zeff)
          effnsq = min(effnsq, 1000.0d0)
          rsqlo  = 2.5d0 * (effnsq / zeff)**2
          if (code - zeff + 1.0d0 > 20.0d0 .and. code - zeff + 1.0d0 < 29.0d0) then
            rsqup = (45.0d0 - (code - zeff + 1.0d0)) / zeff
            rsqlo = 0.0d0
          end if
          block
            character(len=8) :: clabelp_vdw
            clabelp_vdw = transfer(labelp(1), clabelp_vdw)
            if (clabelp_vdw == 'CONTINUU') rsqlo = 0.0d0
          end block
          if (rsqup < rsqlo) rsqup = 2.0d0 * rsqlo
          gammaw_d = 4.5d-9 * (rsqup - rsqlo)**0.4d0
        else
          gammaw_d = 1.d-7 / zeff
        end if
      end if

      ! ---- line type dispatch ----
      itype = 0
      if (code == 1.00)                  itype = -1
      if (code == 1.00 .and. iso1 == 2) itype = -2
      if (code == 2.00)                  itype = -3
      if (code == 2.00 .and. iso1 == 3) itype = -4
      if (code == 2.01)                  itype = -6
      if (code == 2.01 .and. iso1 == 3) itype = -6
      if (auto_flag == 'AUT')            itype =  1
      if (auto_flag == 'PRD')            itype =  3

      block
        character(len=8) :: clabelp
        clabelp = transfer(labelp(1), clabelp)
        if (clabelp == 'CONTINUU') then
          itype = nint(xjp)
          gf    = gf * (xj + xj + 1.0)
        end if
      end block

      ncon = 0
      if (iso1 == 0 .and. iso2 > 0) ncon = iso2

      ! ---- CGF conversion and damping normalisation ----
      if (itype /= 1 .and. itype <= 3) then
        frelin   = 2.99792458d17 / wlvac
        cgf      = 0.026538d0 / 1.77245d0 * dble(gf) / frelin
        frq4pi   = 12.5664d0 * frelin
        if (itype == 2) then
          gammar_d = dble(gr)
        else
          gammar_d = gammar_d / frq4pi
          gammas_d = gammas_d / frq4pi
          gammaw_d = gammaw_d / frq4pi
        end if
      end if

      nbup = abs(nbup)
      nblo = abs(nblo)
      nelionx = 0

      ! ---- append to appropriate buffer ----
      if (itype /= 0) then
        ! complex-profile line -> nlte buffer
        if (nblo + nbup /= 0) then
          do ic = 1, 17
            if (code == codex(ic)) then; nelionx = ic; exit; end if
          end do
          if (nelionx == 0) then
            write(6,'(a,f10.2)') ' WARNING: nlte line has unknown CODE: ', code
            cycle
          end if
        end if
        if (n_nlte == nlte_cap) call grow_nlte(nlte_buf, nlte_cap)
        n_nlte = n_nlte + 1
        if (itype == 1 .or. itype > 3) then
          nlte_buf(n_nlte) = nlte_line_t(wlvac, real(elo_d,4), real(gf,4), &
            nblo, nbup, nelion_i, itype, ncon, nelionx, &
            real(gammar_d,4), real(gammas_d,4), real(gammaw_d,4), nbuff, lim)
        else
          nlte_buf(n_nlte) = nlte_line_t(wlvac, real(elo_d,4), real(cgf,4), &
            nblo, nbup, nelion_i, itype, ncon, nelionx, &
            real(gammar_d,4), real(gammas_d,4), real(gammaw_d,4), nbuff, lim)
        end if

      else if (nblo + nbup /= 0) then
        ! TYPE=0 with departure coefficients -> nlte buffer
        do ic = 1, 17
          if (code == codex(ic)) then; nelionx = ic; exit; end if
        end do
        if (nelionx == 0) then
          write(6,'(a,f10.2)') ' WARNING: nlte line has unknown CODE: ', code
          cycle
        end if
        if (n_nlte == nlte_cap) call grow_nlte(nlte_buf, nlte_cap)
        n_nlte = n_nlte + 1
        nlte_buf(n_nlte) = nlte_line_t(wlvac, real(elo_d,4), real(cgf,4), &
          nblo, nbup, nelion_i, itype, ncon, nelionx, &
          real(gammar_d,4), real(gammas_d,4), real(gammaw_d,4), nbuff, lim)

      else
        ! plain Voigt line -> lte buffer
        if (n_lte == lte_cap) call grow_lte(lte_buf, lte_cap)
        n_lte = n_lte + 1
        lte_buf(n_lte) = lte_line_t(nbuff, real(cgf,4), nelion_i, &
          real(elo_d,4), real(gammar_d,4), real(gammas_d,4), real(gammaw_d,4))
      end if

    end do

    close(11)

    ! ---- trim to exact size and hand back ----
    if (allocated(lte_arr))  deallocate(lte_arr)
    if (allocated(nlte_arr)) deallocate(nlte_arr)
    allocate(lte_arr(n_lte),   source=lte_buf(1:n_lte))
    allocate(nlte_arr(n_nlte), source=nlte_buf(1:n_nlte))
    deallocate(lte_buf, nlte_buf)

  end subroutine read_gfall


  ! ============================================================================
  !  READ_PREDICT — read Kurucz predicted-wavelength binary line list
  !                 (translated from rpredict.for)
  ! ============================================================================
  subroutine read_predict(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                           lte_arr, n_lte)
    character(len=*),             intent(in)    :: filename
    real(8),                      intent(in)    :: wlbeg, wlend, ratiolg
    integer,                      intent(in)    :: ixwlbeg
    type(lte_line_t), allocatable, intent(inout) :: lte_arr(:)
    integer,                      intent(inout) :: n_lte

    integer(4) :: iiiiiii(4)
    integer(4) :: iwl
    integer(2) :: iwords(8)
    equivalence (iiiiiii(1), iwl)
    equivalence (iiiiiii(1), iwords(1))

    integer :: nelionold(1005)
    integer :: nelionolda(209), nelionoldb(286), nelionoldc(95)
    integer :: nelionoldd(95),  nelionolde(95),  nelionoldf(60)
    integer :: nelionoldg(165)
    equivalence (nelionold(  1), nelionolda(1))
    equivalence (nelionold(210), nelionoldb(1))
    equivalence (nelionold(496), nelionoldc(1))
    equivalence (nelionold(591), nelionoldd(1))
    equivalence (nelionold(686), nelionolde(1))
    equivalence (nelionold(781), nelionoldf(1))
    equivalence (nelionold(841), nelionoldg(1))

    data nelionolda/ &
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
    data nelionoldb/ &
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
    data nelionoldc/ &
    181,182,183, 0, 0, 187,188,189, 0, 0, 193,194,195, 0, 0, &
    199,200,201, 0, 0, 205,206,207, 0, 0, 211,212,213, 0, 0, &
    217,218,219, 0, 0, 223,224,225, 0, 0, 229,230,231, 0, 0, &
    235,236,237, 0, 0, 241,242,243, 0, 0, 247,248,249, 0, 0, &
    253,254,255, 0, 0, 259,260,261, 0, 0, 265,266,267, 0, 0, &
    271,272,273, 0, 0, 277,278,279, 0, 0, 283,284,285, 0, 0, &
    289,290,291, 0, 0/
    data nelionoldd/ &
    295,296,297, 0, 0, 301,302,303, 0, 0, 307,308,309, 0, 0, &
    313,314,315, 0, 0, 319,320,321, 0, 0, 325,326,327, 0, 0, &
    331,332,333, 0, 0, 337,338,339, 0, 0, 343,344,345, 0, 0, &
    349,350,351, 0, 0, 355,356,357, 0, 0, 361,362,363, 0, 0, &
    367,368,369, 0, 0, 373,374,375, 0, 0, 379,380,381, 0, 0, &
    385,386,387, 0, 0, 391,392,393, 0, 0, 397,398,399, 0, 0, &
    403,404,405, 0, 0/
    data nelionolde/ &
    409,410,411, 0, 0, 415,416,417, 0, 0, 421,422,423, 0, 0, &
    427,428,429, 0, 0, 433,434,435, 0, 0, 439,440,441, 0, 0, &
    445,446,447, 0, 0, 451,452,453, 0, 0, 457,458,459, 0, 0, &
    463,464,465, 0, 0, 469,470,471, 0, 0, 475,476,477, 0, 0, &
    481,482,483, 0, 0, 487,488,489, 0, 0, 493,494,495, 0, 0, &
    499,500,501, 0, 0, 505,506,507, 0, 0, 511,512,513, 0, 0, &
    517,518,519, 0, 0/
    data nelionoldf/ &
    523,524,525, 0, 0, 529,530,531, 0, 0, 535,536,537, 0, 0, &
    541,542,543, 0, 0, 547,548,549, 0, 0, 553,554,555, 0, 0, &
    559,560,561, 0, 0, 565,566,567, 0, 0, 571,572,573, 0, 0, &
    577,578,579, 0, 0, 583,584,585, 0, 0, 589,590,591, 0, 0/
    data nelionoldg/ &
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

    real(4), save :: tablog(32768)
    real(8) :: wlvac_d, elo_d, gf_d, freq, congf, frq4pi
    real(8) :: gamrf, gamsf, gamwf, ratiolog_r2e6
    integer :: ios, i, n, istart, istop, limitblue, limitred, newlimit
    integer :: lengthfile, nstart, iwl1, nelionnew, nelion_i, nbuff
    integer :: ixwl

    integer, parameter :: CHUNK = 200000
    type(lte_line_t), allocatable :: lte_buf(:)
    integer :: lte_cap

    do i = 1, 32768
      tablog(i) = 10.0**((i - 16384) * 0.001)
    end do

    ratiolog_r2e6 = log(1.0d0 + 1.0d0/2000000.0d0)

    istart = int(log(wlbeg - 1.0d0) / ratiolog_r2e6 + 0.5d0)
    istop  = int(log(wlend + 1.0d0) / ratiolog_r2e6 + 0.5d0)

    open(unit=11, file=filename, status='old', form='unformatted', &
         access='direct', recl=16, iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') 'ERROR: cannot open predict file: ', trim(filename); stop 1
    end if

    read(11, rec=1) iwl1
    if (iwl1 > istop) then
      write(6,'(a)') '  No predicted lines in window (file starts past wlend).'
      close(11)
      if (allocated(lte_arr)) deallocate(lte_arr)
      allocate(lte_arr(0));  n_lte = 0;  return
    end if

    block
      integer(8) :: fsize
      inquire(file=filename, size=fsize)
      lengthfile = int(fsize / 16)
    end block

    read(11, rec=lengthfile) iwl
    if (iwl < istart) then
      write(6,'(a)') '  No predicted lines in window (file ends before wlbeg).'
      close(11)
      if (allocated(lte_arr)) deallocate(lte_arr)
      allocate(lte_arr(0));  n_lte = 0;  return
    end if

    ! binary search for first record >= istart
    limitblue = 1;  limitred = lengthfile
    do while (limitred - limitblue > 1)
      newlimit = (limitred + limitblue) / 2
      read(11, rec=newlimit) iwl
      if (iwl < istart) then
        limitblue = newlimit
      else
        limitred  = newlimit
      end if
    end do
    nstart = limitred

    lte_cap = CHUNK;  allocate(lte_buf(lte_cap))
    n_lte = 0
    n = 0

    do i = nstart, lengthfile
      read(11, rec=i, iostat=ios) iiiiiii
      if (ios < 0) exit
      if (ios > 0) then
        write(6,'(a,i0)') ' ERROR: read_predict: read error at record ', i; exit
      end if
      if (iwl > istop) exit

      nelionnew = abs(iwords(3)) / 10
      if (nelionnew < 1 .or. nelionnew > 1005) cycle
      nelion_i  = nelionold(nelionnew)
      if (nelion_i == 0) cycle

      wlvac_d = exp(iwl * ratiolog_r2e6)
      freq    = 2.99792458d17 / wlvac_d
      gf_d    = tablog(iwords(5))
      congf   = 0.01502d0 * gf_d / freq
      elo_d   = tablog(iwords(4))
      frq4pi  = freq * 12.5664d0
      gamrf   = tablog(iwords(6)) / frq4pi
      gamsf   = tablog(iwords(7)) / frq4pi
      gamwf   = tablog(iwords(8)) / frq4pi

      ixwl  = int(log(wlvac_d) / ratiolg + 0.5d0)
      nbuff = ixwl - ixwlbeg + 1

      if (n_lte == lte_cap) call grow_lte(lte_buf, lte_cap)
      n_lte = n_lte + 1
      lte_buf(n_lte) = lte_line_t(nbuff, real(congf,4), nelion_i, &
        real(elo_d,4), real(gamrf,4), real(gamsf,4), real(gamwf,4))
      n = n + 1
    end do

    close(11)

    if (allocated(lte_arr)) deallocate(lte_arr)
    allocate(lte_arr(n_lte), source=lte_buf(1:n_lte))
    deallocate(lte_buf)

  end subroutine read_predict


  ! ============================================================================
  !  READ_MOLEC — dispatcher: binary sidecar if present, else ASCII
  ! ============================================================================
  subroutine read_molec(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                         lte_arr, n_lte)
    character(len=*),             intent(in)    :: filename
    real(8),                      intent(in)    :: wlbeg, wlend, ratiolg
    integer,                      intent(in)    :: ixwlbeg
    type(lte_line_t), allocatable, intent(inout) :: lte_arr(:)
    integer,                      intent(inout) :: n_lte

    character(len=512) :: binfile
    logical :: have_bin
    integer :: dot

    ! Build sidecar path: replace extension with .bin
    binfile = trim(filename)
    dot     = index(binfile, '.', back=.true.)
    if (dot > 0) then
      binfile = binfile(:dot-1) // '.bin'
    else
      binfile = trim(binfile) // '.bin'
    end if
    inquire(file=binfile, exist=have_bin)

    if (have_bin) then
      call read_molec_bin(filename, binfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                          lte_arr, n_lte)
    else
      call read_molec_ascii(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                            lte_arr, n_lte)
    end if

  end subroutine read_molec


  ! ============================================================================
  !  READ_MOLEC_BIN — binary sidecar path with O(log N) window seek
  ! ============================================================================
  subroutine read_molec_bin(ascfile, binfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                             lte_arr, n_lte)
    character(len=*),             intent(in)    :: ascfile, binfile
    real(8),                      intent(in)    :: wlbeg, wlend, ratiolg
    integer,                      intent(in)    :: ixwlbeg
    type(lte_line_t), allocatable, intent(inout) :: lte_arr(:)
    integer,                      intent(inout) :: n_lte

    integer, parameter :: MAGIC    = int(z'4D4F4C32')
    integer, parameter :: RECBYTES = 32

    integer(4) :: hdr(8), buf(8)
    real(4)    :: wlvac, gflog, e_r4, ep_r4
    integer(4) :: icode, iso, loggr, labelp_x
    integer    :: nrec, ios, ixwl, nbuff, lo, hi, mid, istart
    real(8)    :: wlvac_d, elo, gf, freq, congf, frq4pi
    real(8)    :: gammar, gammas, gammaw, gamrf, gamsf, gamwf
    real(4)    :: fudge, x1, x2
    integer    :: nelion, iso1, iso2

    integer, parameter :: CHUNK = 200000
    type(lte_line_t), allocatable :: lte_buf(:)
    ! n_lte is a running total across molecule files; preserve existing contents
    type(lte_line_t), allocatable :: lte_old(:)
    integer :: n_old, n_new, lte_cap

    open(unit=21, file=binfile, status='old', form='unformatted', &
         access='direct', recl=RECBYTES, iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') '  WARNING: cannot open binary sidecar, falling back: ', &
        trim(binfile)
      call read_molec_ascii(ascfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                            lte_arr, n_lte)
      return
    end if

    read(21, rec=1) hdr
    if (hdr(1) /= MAGIC) then
      write(6,'(a)') '  WARNING: binary sidecar wrong magic, falling back'
      close(21)
      call read_molec_ascii(ascfile, wlbeg, wlend, ratiolg, ixwlbeg, &
                            lte_arr, n_lte)
      return
    end if
    nrec = hdr(2)

    ! Binary search: first record with wlvac >= wlbeg - 1.0
    lo = 1;  hi = nrec
    do while (lo < hi)
      mid = (lo + hi) / 2
      read(21, rec=mid+1, iostat=ios) buf
      if (ios /= 0) then
        write(6,'(a,i0)') ' ERROR: read_molec_bin seek error at record ', mid+1
        close(21);  return
      end if
      wlvac = transfer(buf(1), wlvac)
      if (wlvac < real(wlbeg - 1.0d0, 4)) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    istart = lo

    lte_cap = CHUNK;  allocate(lte_buf(lte_cap))
    n_new = 0

    do mid = istart, nrec
      read(21, rec=mid+1, iostat=ios) buf
      if (ios < 0) exit
      if (ios > 0) then
        write(6,'(a,i0)') ' ERROR: read_molec_bin read error at record ', mid+1
        exit
      end if

      wlvac    = transfer(buf(1), wlvac)
      gflog    = transfer(buf(2), gflog)
      e_r4     = transfer(buf(3), e_r4)
      ep_r4    = transfer(buf(4), ep_r4)
      icode    = buf(5)
      iso      = buf(6)
      loggr    = buf(7)
      labelp_x = buf(8)

      wlvac_d = dble(wlvac)
      if (wlvac_d > wlend + 1.0d0) exit

      call molec_dispatch(iso, icode, nelion, iso1, iso2, x1, x2, fudge)
      if (nelion == 0) cycle

      gf  = exp((dble(gflog) + x1 + x2 + fudge) * 2.30258509299405d0)
      elo = min(abs(dble(e_r4)), abs(dble(ep_r4)))

      ixwl   = int(log(wlvac_d) / ratiolg + 0.5d0)
      nbuff  = ixwl - ixwlbeg + 1
      freq   = 2.99792458d17 / wlvac_d
      congf  = 0.026538d0 / 1.77245d0 * gf / freq
      frq4pi = freq * 12.5664d0

      if (loggr == 0) then
        gammar = 2.223d13 / wlvac_d**2
      else
        gammar = 10.0d0**(loggr * 0.01d0)
      end if

      if (labelp_x == 1) then
        gammas = 3.0d-8;  gammaw = 1.0d-8
      else
        gammas = 3.0d-5;  gammaw = 1.0d-7
      end if

      gamrf = gammar / frq4pi
      gamsf = gammas / frq4pi
      gamwf = gammaw / frq4pi

      if (n_new == lte_cap) call grow_lte(lte_buf, lte_cap)
      n_new = n_new + 1
      lte_buf(n_new) = lte_line_t(nbuff, real(congf,4), nelion, &
        real(elo,4), real(gamrf,4), real(gamsf,4), real(gamwf,4))
    end do

    close(21)

    ! Append to existing lte_arr (which may already have lines from prior mol files)
    n_old = n_lte
    if (n_old > 0 .and. allocated(lte_arr)) then
      allocate(lte_old(n_old), source=lte_arr(1:n_old))
    end if
    if (allocated(lte_arr)) deallocate(lte_arr)
    n_lte = n_old + n_new
    allocate(lte_arr(n_lte))
    if (n_old > 0) lte_arr(1:n_old) = lte_old(1:n_old)
    if (n_new > 0) lte_arr(n_old+1:n_lte) = lte_buf(1:n_new)
    if (allocated(lte_old)) deallocate(lte_old)
    deallocate(lte_buf)

  end subroutine read_molec_bin


  ! ============================================================================
  !  READ_MOLEC_ASCII — sequential ASCII scan with window guards
  ! ============================================================================
  subroutine read_molec_ascii(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                               lte_arr, n_lte)
    character(len=*),             intent(in)    :: filename
    real(8),                      intent(in)    :: wlbeg, wlend, ratiolg
    integer,                      intent(in)    :: ixwlbeg
    type(lte_line_t), allocatable, intent(inout) :: lte_arr(:)
    integer,                      intent(inout) :: n_lte

    real(8)          :: wl, e, ep, wlvac, elo, gf, freq, congf, frq4pi
    real(8)          :: gammar, gammas, gammaw, gamrf, gamsf, gamwf
    real(4)          :: gflog, xj, xjp, fudge, x1, x2
    integer          :: loggr, icode, iso, nelion, iso1, iso2
    integer          :: ios, ixwl, nbuff
    character(len=8) :: clabel, clabelp

    integer, parameter :: CHUNK = 200000
    type(lte_line_t), allocatable :: lte_buf(:)
    type(lte_line_t), allocatable :: lte_old(:)
    integer :: n_old, n_new, lte_cap

    open(unit=11, file=filename, status='old', action='read', &
         form='formatted', iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') 'ERROR: cannot open molecule file: ', trim(filename); stop 1
    end if

    lte_cap = CHUNK;  allocate(lte_buf(lte_cap))
    n_new = 0

    do
      read(11, '(F10.4,F7.3,F5.1,F10.3,F5.1,F11.3,I4,A8,A8,I2,I4)', &
           iostat=ios) wl, gflog, xj, e, xjp, ep, icode, clabel, clabelp, iso, loggr
      if (ios < 0) exit
      if (ios > 0) then
        write(6,'(a,i0)') ' ERROR: read_molec_ascii: read error near line ', n_new; exit
      end if

      if (abs(wl) > wlend + 2.0d0) exit

      wlvac = 1.0d7 / abs(abs(ep) - abs(e))
      if (wlvac < wlbeg - 1.0d0) cycle
      if (wlvac > wlend + 1.0d0) cycle

      call molec_dispatch(iso, icode, nelion, iso1, iso2, x1, x2, fudge)
      if (nelion == 0) cycle

      gf  = exp((gflog + x1 + x2 + fudge) * 2.30258509299405d0)
      elo = min(abs(e), abs(ep))

      ixwl   = int(log(wlvac) / ratiolg + 0.5d0)
      nbuff  = ixwl - ixwlbeg + 1
      freq   = 2.99792458d17 / wlvac
      congf  = 0.026538d0 / 1.77245d0 * gf / freq
      frq4pi = freq * 12.5664d0

      if (loggr == 0) then
        gammar = 2.223d13 / wlvac**2
      else
        gammar = 10.0d0**(loggr * 0.01d0)
      end if

      if (clabelp(1:1) == 'X') then
        gammas = 3.0d-8;  gammaw = 1.0d-8
      else
        gammas = 3.0d-5;  gammaw = 1.0d-7
      end if

      gamrf = gammar / frq4pi
      gamsf = gammas / frq4pi
      gamwf = gammaw / frq4pi

      if (n_new == lte_cap) call grow_lte(lte_buf, lte_cap)
      n_new = n_new + 1
      lte_buf(n_new) = lte_line_t(nbuff, real(congf,4), nelion, &
        real(elo,4), real(gamrf,4), real(gamsf,4), real(gamwf,4))
    end do

    close(11)

    n_old = n_lte
    if (n_old > 0 .and. allocated(lte_arr)) then
      allocate(lte_old(n_old), source=lte_arr(1:n_old))
    end if
    if (allocated(lte_arr)) deallocate(lte_arr)
    n_lte = n_old + n_new
    allocate(lte_arr(n_lte))
    if (n_old > 0) lte_arr(1:n_old) = lte_old(1:n_old)
    if (n_new > 0) lte_arr(n_old+1:n_lte) = lte_buf(1:n_new)
    if (allocated(lte_old)) deallocate(lte_old)
    deallocate(lte_buf)

  end subroutine read_molec_ascii


  ! ============================================================================
  !  READ_H2O — read Partridge-Schwenke H2O binary line list
  !              (translated from rh2ofast.for)
  ! ============================================================================
  subroutine read_h2o(filename, wlbeg, wlend, ratiolg, ixwlbeg, &
                       lte_arr, n_lte)
    character(len=*),             intent(in)    :: filename
    real(8),                      intent(in)    :: wlbeg, wlend, ratiolg
    integer,                      intent(in)    :: ixwlbeg
    type(lte_line_t), allocatable, intent(inout) :: lte_arr(:)
    integer,                      intent(inout) :: n_lte

    integer(4) :: irec(2)
    integer(4) :: iwl
    integer(2) :: ielo_i
    equivalence (irec(1), iwl)
    equivalence (irec(2), ielo_i)

    real(4), save :: tablog(32768)
    real(4), parameter :: xiso(4)  = [ 0.9976,  0.0004,  0.0020, 0.00001]
    real(4), parameter :: x2iso(4) = [-0.001,  -3.398,  -2.690, -5.000]

    real(8) :: wlvac, wlvac1, freq, congf, ratiolog_r2e6, frq4pi
    real(8) :: gammar, gamrf, gamsf, gamwf
    integer :: ios, iline, istart, i, iso, ixwl, nbuff
    integer :: limitblue, limitred, newlimit, lengthfile
    integer(2) :: igflog_loc, ielo_loc

    integer, parameter :: CHUNK = 500000
    type(lte_line_t), allocatable :: lte_buf(:)
    integer :: n_new, lte_cap

    do i = 1, 32768
      tablog(i) = 10.0**((i - 16384) * 0.001)
    end do

    ratiolog_r2e6 = log(1.0d0 + 1.0d0/2000000.0d0)

    open(unit=11, file=filename, status='old', form='unformatted', &
         access='direct', recl=8, iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') 'ERROR: cannot open H2O file: ', trim(filename); stop 1
    end if

    read(11, rec=1) irec
    wlvac = exp(iwl * ratiolog_r2e6)
    if (wlvac > wlend + 1.0d0) then
      close(11)
      if (.not. allocated(lte_arr)) allocate(lte_arr(0))
      return
    end if

    block
      integer(8) :: fsize
      inquire(file=filename, size=fsize)
      lengthfile = int(fsize / 8)
    end block

    read(11, rec=lengthfile) irec
    wlvac1 = exp(iwl * ratiolog_r2e6)
    if (wlbeg - 1.0d0 > wlvac1) then
      close(11)
      if (.not. allocated(lte_arr)) allocate(lte_arr(0))
      return
    end if

    limitblue = 1;  limitred = lengthfile
    do while (limitred - limitblue > 1)
      newlimit = (limitred + limitblue) / 2
      read(11, rec=newlimit) irec
      wlvac = exp(iwl * ratiolog_r2e6)
      if (wlvac < wlbeg - 1.0d0) then
        limitblue = newlimit
      else
        limitred  = newlimit
      end if
    end do
    istart = newlimit

    lte_cap = CHUNK;  allocate(lte_buf(lte_cap))
    n_new  = 0
    iline  = 0

    do i = istart, lengthfile
      read(11, rec=i, iostat=ios) irec
      if (ios < 0) exit
      if (ios > 0) then
        write(6,'(a,i0)') ' ERROR: read_h2o: read error at record ', i; exit
      end if

      ielo_loc   = int(ibits(irec(2),  0, 16), 2)
      igflog_loc = int(ibits(irec(2), 16, 16), 2)

      wlvac = exp(iwl * ratiolog_r2e6)
      if (wlvac > wlend + 1.0d0) exit

      freq  = 2.99792458d17 / wlvac
      ixwl  = int(log(wlvac) / ratiolg + 0.5d0)
      nbuff = ixwl - ixwlbeg + 1

      if      (ielo_loc > 0 .and. igflog_loc > 0)  then; iso = 1
      else if (ielo_loc > 0 .and. igflog_loc <= 0) then; iso = 2
      else if (ielo_loc <= 0 .and. igflog_loc > 0) then; iso = 3
      else;                                               iso = 4
      end if

      congf  = 0.01502d0 * tablog(abs(igflog_loc)) / freq * xiso(iso)
      frq4pi = freq * 12.5664d0
      gammar = 2.223d13 / wlvac**2 * 0.001d0
      gamrf  = gammar / frq4pi
      gamsf  = tablog(1)    / frq4pi
      gamwf  = tablog(9384) / frq4pi

      if (n_new == lte_cap) call grow_lte(lte_buf, lte_cap)
      n_new = n_new + 1
      lte_buf(n_new) = lte_line_t(nbuff, real(congf,4), 534, &
        real(abs(ielo_loc),4), real(gamrf,4), real(gamsf,4), real(gamwf,4))
      iline = iline + 1
    end do

    close(11)

    if (allocated(lte_arr)) deallocate(lte_arr)
    allocate(lte_arr(n_new), source=lte_buf(1:n_new))
    n_lte = n_new
    deallocate(lte_buf)

  end subroutine read_h2o


  ! ============================================================================
  !  MOLEC_DISPATCH — map ISO isotope index to NELION and abundance corrections
  ! ============================================================================
  subroutine molec_dispatch(iso, icode, nelion, iso1, iso2, x1, x2, fudge)
    integer, intent(in)  :: iso, icode
    integer, intent(out) :: nelion, iso1, iso2
    real(4), intent(out) :: x1, x2, fudge

    nelion = 0; iso1 = 0; iso2 = 0; x1 = 0.0; x2 = 0.0; fudge = 0.0

    select case (iso)
    case (1);   nelion=240; iso1=1;  iso2=1;  x1= 0.0;    x2=-5.0
    case (2);   nelion=240; iso1=1;  iso2=2;  x1= 0.0;    x2=-4.469
    case (12)
      select case (icode)
      case (606); nelion=264; iso1=12; iso2=12; x1=-.005; x2=-.005
      case (608); nelion=276; iso1=12; iso2=16; x1=-.005; x2=-.001
      case (106); nelion=246; iso1=1;  iso2=12; x1= 0.0;  x2=-.005
      case default
                  nelion=270; iso1=12; iso2=14; x1=-.005; x2=-.002
      end select
    case (13)
      select case (icode)
      case (606); nelion=264; iso1=12; iso2=13; x1=-.005;  x2=-1.955
      case (608); nelion=276; iso1=13; iso2=16; x1=-1.955; x2=-.001
      case (106); nelion=246; iso1=1;  iso2=13; x1= 0.0;   x2=-1.955
      case default
                  nelion=270; iso1=13; iso2=14; x1=-1.955; x2=-.002
      end select
    case (14);  nelion=252; iso1=1;  iso2=14; x1= 0.0;   x2=-.002
    case (15)
      if (icode == 607) then
                  nelion=270; iso1=12; iso2=15; x1=-.005; x2=-2.444
      else;       nelion=252; iso1=1;  iso2=15; x1= 0.0;  x2=-2.444
      end if
    case (16)
      if (icode == 813) then
                  nelion=324; iso1=27; iso2=16; x1= 0.0;  x2=-.001
      else;       nelion=258; iso1=1;  iso2=16; x1= 0.0;  x2=-.001
      end if
    case (17)
      if (icode == 813) then
                  nelion=324; iso1=27; iso2=17; x1= 0.0;  x2=-3.398
      else;       nelion=276; iso1=12; iso2=17; x1=-.005; x2=-3.398
      end if
    case (18)
      select case (icode)
      case (814); nelion=330; iso1=28; iso2=18; x1=-.035; x2=-2.690
      case (608); nelion=276; iso1=12; iso2=18; x1=-.005; x2=-2.690
      case (813); nelion=324; iso1=27; iso2=18; x1= 0.0;  x2=-2.690
      case default
                  nelion=258; iso1=1;  iso2=18; x1= 0.0;  x2=-2.690
      end select
    case (23);  nelion=492; iso1=1;  iso2=23; x1= 0.0;   x2= 0.0
    case (24)
      if (icode == 812) then
                  nelion=318; iso1=16; iso2=24; x1=-.001; x2=-.102
      else;       nelion=300; iso1=1;  iso2=24; x1= 0.0;  x2=-.105
      end if
    case (25)
      if (icode == 812) then
                  nelion=318; iso1=16; iso2=25; x1=-.001; x2=-1.000
      else;       nelion=300; iso1=1;  iso2=25; x1= 0.0;  x2=-.996
      end if
    case (26)
      if (icode == 812) then
                  nelion=318; iso1=16; iso2=26; x1=-.001; x2=-.958
      else;       nelion=300; iso1=1;  iso2=26; x1= 0.0;  x2=-.947
      end if
    case (28)
      if (icode == 814) then
                  nelion=330; iso1=28; iso2=16; x1=-.035; x2=-.001
      else;       nelion=312; iso1=1;  iso2=28; x1= 0.0;  x2=-.035
      end if
    case (29)
      if (icode == 814) then
                  nelion=330; iso1=29; iso2=16; x1=-1.328; x2=-.001
      else;       nelion=312; iso1=1;  iso2=29; x1= 0.0;   x2=-1.331
      end if
    case (30)
      if (icode == 814) then
                  nelion=330; iso1=30; iso2=16; x1=-1.510; x2=-.001
      else;       nelion=312; iso1=1;  iso2=30; x1= 0.0;   x2=-1.516
      end if
    case (33);  nelion=264; iso1=13; iso2=13; x1=-1.955; x2=-1.955
    case (39);  nelion=498; iso1=39; iso2=1;  x1=-.030;  x2= 0.0
    case (40)
      if (icode == 820) then
                  nelion=354; iso1=40; iso2=16; x1=-.013; x2=-.001
      else;       nelion=342; iso1=40; iso2=1;  x1=-.013; x2= 0.0
      end if
    case (41);  nelion=498; iso1=41; iso2=1;  x1=-1.172; x2= 0.0
    case (42);  nelion=342; iso1=42; iso2=1;  x1=-2.189; x2= 0.0
    case (43);  nelion=342; iso1=43; iso2=1;  x1=-2.870; x2= 0.0
    case (44);  nelion=342; iso1=44; iso2=1;  x1=-1.681; x2= 0.0
    case (46)
      if (icode == 120) then
                  nelion=342; iso1=46; iso2=1;  x1=-4.398; x2= 0.0
      else;       nelion=366; iso1=16; iso2=46; x1= 0.0;   x2=-1.101
      end if
    case (47);  nelion=366; iso1=16; iso2=47; x1= 0.0;   x2=-1.138
    case (48)
      if (icode == 120) then
                  nelion=342; iso1=48; iso2=1;  x1=-2.728; x2= 0.0
      else;       nelion=366; iso1=16; iso2=48; x1= 0.0;   x2=-.131
      end if
    case (49);  nelion=366; iso1=16; iso2=49; x1= 0.0;   x2=-1.259
    case (50)
      if (icode == 124) then
                  nelion=432; iso1=50; iso2=1;  x1=-1.362; x2= 0.0
      else;       nelion=366; iso1=16; iso2=50; x1= 0.0;   x2=-1.272
      end if
    case (51);  nelion=372; iso1=16; iso2=51; x1= 0.0;   x2=-.001
    case (52);  nelion=432; iso1=52; iso2=1;  x1=-.077;  x2= 0.0
    case (53);  nelion=432; iso1=53; iso2=1;  x1=-1.022; x2= 0.0
    case (54)
      if (icode == 124) then
                  nelion=432; iso1=54; iso2=1;  x1=-1.626; x2= 0.0
      else;       nelion=444; iso1=54; iso2=1;  x1=-1.237; x2= 0.0
      end if
    case (56);  nelion=444; iso1=56; iso2=1;  x1=-.038;  x2= 0.0
    case (57);  nelion=444; iso1=57; iso2=1;  x1=-1.658; x2= 0.0
    case (58);  nelion=444; iso1=58; iso2=1;  x1=-2.553; x2= 0.0
    case default
      nelion = 0
    end select
  end subroutine molec_dispatch


  ! ============================================================================
  !  GROW_LTE / GROW_NLTE — double buffer capacity
  ! ============================================================================
  subroutine grow_lte(buf, cap)
    type(lte_line_t), allocatable, intent(inout) :: buf(:)
    integer,                       intent(inout) :: cap
    type(lte_line_t), allocatable :: tmp(:)
    allocate(tmp(cap*2))
    tmp(1:cap) = buf(1:cap)
    call move_alloc(tmp, buf)
    cap = cap * 2
  end subroutine grow_lte

  subroutine grow_nlte(buf, cap)
    type(nlte_line_t), allocatable, intent(inout) :: buf(:)
    integer,                        intent(inout) :: cap
    type(nlte_line_t), allocatable :: tmp(:)
    allocate(tmp(cap*2))
    tmp(1:cap) = buf(1:cap)
    call move_alloc(tmp, buf)
    cap = cap * 2
  end subroutine grow_nlte


  ! ============================================================================
  !  IONPOT_INDEX — upper-level effective principal quantum number
  ! ============================================================================
  subroutine ionpot_index(nelem, icharge, eup, effnsq, zeff)
    integer, intent(in)    :: nelem, icharge
    real(8), intent(in)    :: eup, zeff
    real(8), intent(inout) :: effnsq   ! in: safety default; out: computed value
    integer :: idx
    real(8) :: deleup

    if (nelem <= 30) then
      idx = nelem*(nelem+1)/2 + icharge
    else
      idx = nelem*5 + 341 + icharge
    end if
    deleup = potion(idx) - eup
    if (deleup > 0.0d0) then
      effnsq = 109737.31d0 * zeff**2 / deleup
    end if
    ! else: retain the caller-supplied safety default (25.0)
  end subroutine ionpot_index


  ! ============================================================================
  !  IONPOT_INDEX_LO — lower-level effective principal quantum number
  ! ============================================================================
  subroutine ionpot_index_lo(nelem, icharge, elo, effnsq, zeff)
    integer, intent(in)  :: nelem, icharge
    real(8), intent(in)  :: elo, zeff
    real(8), intent(out) :: effnsq
    integer :: idx
    real(8) :: delelo

    if (nelem <= 30) then
      idx = nelem*(nelem+1)/2 + icharge
    else
      idx = nelem*5 + 341 + icharge
    end if
    delelo = potion(idx) - elo
    effnsq = 109737.31d0 * zeff**2 / delelo
    effnsq = min(effnsq, 1000.0d0)
  end subroutine ionpot_index_lo



  ! ============================================================================
  !  BASENAME — return the filename portion of a path (after last '/')
  ! ============================================================================
  pure function basename(path) result(name)
    character(len=*), intent(in) :: path
    character(len=512)           :: name
    integer :: i
    i = index(path, '/', back=.true.)
    if (i > 0) then
      name = path(i+1:)
    else
      name = path
    end if
  end function basename


  ! ============================================================================
  !  IS_TIO_FILE — heuristic: does the basename begin with "tio" (any case)?
  !
  !  Used to gate TiO line lists out of hot-star runs.  Kurucz molecular line
  !  list filenames conventionally lead with the species symbol (tioXXX.asc,
  !  tiopred.bin, etc.), so a prefix match is both specific enough to avoid
  !  false positives (e.g. "ratio", "station") and general enough to cover
  !  the filenames Kurucz distributes.
  ! ============================================================================
  pure function is_tio_file(path) result(isTiO)
    character(len=*), intent(in) :: path
    logical                      :: isTiO
    character(len=512) :: bname
    character(len=3)   :: prefix
    integer :: k, c

    bname = basename(path)
    prefix = bname(1:3)
    ! Lowercase in place (ASCII only; Kurucz filenames are ASCII).
    do k = 1, 3
      c = iachar(prefix(k:k))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        prefix(k:k) = achar(c + 32)
      end if
    end do
    isTiO = (prefix == 'tio')
  end function is_tio_file


  ! ============================================================================
  !  IONPOTS — fill the ionisation potential table (cm^-1)
  ! ============================================================================
  subroutine ionpots()
    potion = 0.0d0
    ! H, He
    potion(  1) = 109678.772d0
    potion(  3) = 198310.666d0; potion(  4) = 438908.879d0
    ! Li
    potion(  6) =  43487.114d0; potion(  7) = 610078.526d0
    potion(  8) = 987661.014d0
    ! Be
    potion( 10) =  75192.640d0; potion( 11) = 146882.86d0
    potion( 12) = 1241256.600d0; potion( 13) = 1756018.822d0
    ! B
    potion( 15) =  66928.040d0; potion( 16) = 202887.40d0
    potion( 17) = 305930.80d0;  potion( 18) = 2091972.d0
    potion( 19) = 2744107.936d0
    ! C
    potion( 21) =  90820.42d0;  potion( 22) = 196674.d0
    potion( 23) = 386241.0d0;   potion( 24) = 520175.8d0
    potion( 25) = 3162423.30d0; potion( 26) = 3952061.670d0
    ! N
    potion( 28) = 117225.70d0;  potion( 29) = 238750.20d0
    potion( 30) = 382672.d0;    potion( 31) = 624866.d0
    potion( 32) = 789537.d0;    potion( 33) = 4452723.30d0
    potion( 34) = 5380089.80d0
    ! O
    potion( 36) = 109837.02d0;  potion( 37) = 283270.90d0
    potion( 38) = 443085.0d0;   potion( 39) = 624382.0d0
    potion( 40) = 918657.d0;    potion( 41) = 1114004.d0
    potion( 42) = 5963073.00d0; potion( 43) = 7028394.70d0
    ! F
    potion( 45) = 140524.50d0;  potion( 46) = 282058.6d0
    potion( 47) = 505774.0d0;   potion( 48) = 703110.d0
    potion( 49) = 921480.d0;    potion( 50) = 1267606.0d0
    potion( 51) = 1493632.d0;   potion( 52) = 7693706.60d0
    potion( 53) = 8897242.50d0
    ! Ne
    potion( 55) = 173929.750d0; potion( 56) = 330388.60d0
    potion( 57) = 511544.d0;    potion( 58) = 783890.d0
    potion( 59) = 1018250.d0;   potion( 60) = 1273820.d0
    potion( 61) = 1671750.d0;   potion( 62) = 1928447.d0
    potion( 63) = 9644840.7d0;  potion( 64) = 10986877.20d0
    ! Na
    potion( 66) =  41449.451d0; potion( 67) = 381390.2d0
    potion( 68) = 577654.d0;    potion( 69) = 797970.d0
    potion( 70) = 1116300.d0;   potion( 71) = 1389100.d0
    potion( 72) = 1681700.d0;   potion( 73) = 2130850.d0
    potion( 74) = 2418500.d0;   potion( 75) = 11817106.70d0
    potion( 76) = 13297680.0d0
    ! Mg
    potion( 78) =  61671.050d0; potion( 79) = 121267.61d0
    potion( 80) = 646402.d0;    potion( 81) = 881285.d0
    potion( 82) = 1139900.d0;   potion( 83) = 1506300.d0
    potion( 84) = 1814900.d0;   potion( 85) = 2144820.d0
    potion( 86) = 2645400.d0;   potion( 87) = 2964000.d0
    potion( 88) = 14209914.7d0; potion( 89) = 15829950.d0
    ! Al
    potion( 91) =  48278.48d0;  potion( 92) = 151862.50d0
    potion( 93) = 229445.70d0;  potion( 94) = 967804.d0
    potion( 95) = 1240684.d0;   potion( 96) = 1536400.d0
    potion( 97) = 1949900.d0;   potion( 98) = 2295800.d0
    potion( 99) = 2663300.d0;   potion(100) = 3215300.d0
    potion(101) = 3565010.d0;   potion(102) = 16824539.3d0
    potion(103) = 18584143.0d0
    ! Si
    potion(105) =  65747.76d0;  potion(106) = 131838.14d0
    potion(107) = 270139.30d0;  potion(108) = 364093.10d0
    potion(109) = 1345070.d0;   potion(110) = 1655590.d0
    potion(111) = 1986700.d0;   potion(112) = 2449200.d0
    potion(113) = 2831800.d0;   potion(114) = 3237400.d0
    potion(115) = 3840600.d0;   potion(116) = 4221630.d0
    potion(117) = 19661038.9d0; potion(118) = 21560631.0d0
    ! P
    potion(120) =  84580.83d0;  potion(121) = 159451.70d0
    potion(122) = 243600.70d0;  potion(123) = 414922.8d0
    potion(124) = 524462.9d0;   potion(125) = 1777890.d0
    potion(126) = 2125800.d0;   potion(127) = 2497100.d0
    potion(128) = 3002900.d0;   potion(129) = 3423000.d0
    potion(130) = 3867000.d0;   potion(131) = 4521700.d0
    potion(132) = 4934020.d0;   potion(133) = 22719901.6d0
    potion(134) = 24759942.d0
    ! S
    potion(136) =  83559.1d0;   potion(137) = 188232.7d0
    potion(138) = 281100.d0;    potion(139) = 380870.d0
    potion(140) = 585514.d0;    potion(141) = 710195.d0
    potion(142) = 2266050.d0;   potion(143) = 2651900.d0
    potion(144) = 3063600.d0;   potion(145) = 3611300.d0
    potion(146) = 4069500.d0;   potion(147) = 4552200.d0
    potion(148) = 5258400.d0;   potion(149) = 5702290.d0
    potion(150) = 26001545.1d0; potion(151) = 28182526.d0
    ! Cl
    potion(153) = 104591.00d0;  potion(154) = 192070.0d0
    potion(155) = 321000.d0;    potion(156) = 429400.d0
    potion(157) = 545800.d0;    potion(158) = 781900.d0
    potion(159) = 921096.d0;    potion(160) = 2809280.d0
    potion(161) = 3233080.d0;   potion(162) = 3683000.d0
    potion(163) = 4274000.d0;   potion(164) = 4771400.d0
    potion(165) = 5293400.d0;   potion(166) = 6051000.d0
    potion(167) = 6526620.d0;   potion(168) = 29506532.5d0
    potion(169) = 31828983.d0
    ! Ar
    potion(171) = 127109.842d0; potion(172) = 222848.3d0
    potion(173) = 328550.d0;    potion(174) = 480600.d0
    potion(175) = 603700.d0;    potion(176) = 736300.d0
    potion(177) = 1003400.d0;   potion(178) = 1157056.d0
    potion(179) = 3408500.d0;   potion(180) = 3869500.d0
    potion(181) = 4359000.d0;   potion(182) = 4992000.d0
    potion(183) = 5528700.d0;   potion(184) = 6090500.d0
    potion(185) = 6899800.d0;   potion(186) = 7407190.d0
    potion(187) = 33235410.d0;  potion(188) = 35699895.d0
    ! K
    potion(190) =  35009.814d0; potion(191) = 255072.8d0
    potion(192) = 369427.d0;    potion(193) = 491330.d0
    potion(194) = 666700.d0;    potion(195) = 802000.d0
    potion(196) = 948200.d0;    potion(197) = 1249100.d0
    potion(198) = 1418063.d0;   potion(199) = 4062400.d0
    potion(200) = 4562000.d0;   potion(201) = 5090000.d0
    potion(202) = 5764000.d0;   potion(203) = 6342000.d0
    potion(204) = 6943800.d0;   potion(205) = 7805000.d0
    potion(206) = 8344140.d0;   potion(207) = 37189176.0d0
    potion(208) = 39795784.d0
    ! Ca
    potion(210) =  49305.924d0; potion(211) =  95751.870d0
    potion(212) = 410642.3d0;   potion(213) = 542595.d0
    potion(214) = 680200.d0;    potion(215) = 877400.d0
    potion(216) = 1026000.d0;   potion(217) = 1187600.d0
    potion(218) = 1520600.d0;   potion(219) = 1704050.d0
    potion(220) = 4771600.d0;   potion(221) = 5309000.d0
    potion(222) = 5877000.d0;   potion(223) = 6591000.d0
    potion(224) = 7210000.d0;   potion(225) = 7853000.d0
    potion(226) = 8766000.d0;   potion(227) = 9337690.d0
    potion(228) = 41367028.d0;  potion(229) = 44117409.d0
    ! Sc
    potion(231) =  52922.00d0;  potion(232) = 103237.1d0
    potion(233) = 199677.37d0;  potion(234) = 592732.d0
    potion(235) = 741600.d0;    potion(236) = 892700.d0
    potion(237) = 1113000.d0;   potion(238) = 1275000.d0
    potion(239) = 1452000.d0;   potion(240) = 1816200.d0
    potion(241) = 2014760.d0;   potion(242) = 5543900.d0
    potion(243) = 6111000.d0;   potion(244) = 6720000.d0
    potion(245) = 7473000.d0;   potion(246) = 8135000.d0
    potion(247) = 8820000.d0;   potion(248) = 9784000.d0
    potion(249) = 10388070.d0;  potion(250) = 45771185.d0
    potion(251) = 48665510.d0
    ! Ti
    potion(253) =  55072.50d0;  potion(254) = 109494.d0
    potion(255) = 221735.6d0;   potion(256) = 348973.3d0
    potion(257) = 800900.d0;    potion(258) = 964100.d0
    potion(259) = 1134700.d0;   potion(260) = 1375000.d0
    potion(261) = 1549000.d0;   potion(262) = 1741500.d0
    potion(263) = 2137900.d0;   potion(264) = 2351110.d0
    potion(265) = 6353000.d0;   potion(266) = 6969000.d0
    potion(267) = 7618000.d0;   potion(268) = 8408000.d0
    potion(269) = 9116000.d0;   potion(270) = 9842000.d0
    potion(271) = 10859000.d0;  potion(272) = 11495470.d0
    potion(273) = 50401766.d0;  potion(274) = 53440740.d0
    ! V
    potion(276) =  54411.67d0;  potion(277) = 117900.d0
    potion(278) = 236410.d0;    potion(279) = 376730.d0
    potion(280) = 526532.0d0;   potion(281) = 1033400.d0
    potion(282) = 1215700.d0;   potion(283) = 1399800.d0
    potion(284) = 1661000.d0;   potion(285) = 1859000.d0
    potion(286) = 2055000.d0;   potion(287) = 2488200.d0
    potion(288) = 2712230.d0;   potion(289) = 7227000.d0
    potion(290) = 7882000.d0;   potion(291) = 8573000.d0
    potion(292) = 9398000.d0;   potion(293) = 10153000.d0
    potion(294) = 10922000.d0;  potion(295) = 11991000.d0
    potion(296) = 12660130.d0;  potion(297) = 55259549.d0
    potion(298) = 58443920.d0
    ! Cr
    potion(300) =  54575.6d0;   potion(301) = 132971.02d0
    potion(302) = 249700.d0;    potion(303) = 396500.d0
    potion(304) = 560200.d0;    potion(305) = 731020.d0
    potion(306) = 1292800.d0;   potion(307) = 1490200.d0
    potion(308) = 1690100.d0;   potion(309) = 1972000.d0
    potion(310) = 2184000.d0;   potion(311) = 2393000.d0
    potion(312) = 2860500.d0;   potion(313) = 3098480.d0
    potion(314) = 8159000.d0;   potion(315) = 8850000.d0
    potion(316) = 9582000.d0;   potion(317) = 10443000.d0
    potion(318) = 11247000.d0;  potion(319) = 12059000.d0
    potion(320) = 13180000.d0;  potion(321) = 13882280.d0
    potion(322) = 60345293.d0;  potion(323) = 63675850.d0
    ! Mn
    potion(325) =  59959.4d0;   potion(326) = 126145.00d0
    potion(327) = 271550.d0;    potion(328) = 413000.d0
    potion(329) = 584000.d0;    potion(330) = 771100.d0
    potion(331) = 961440.d0;    potion(332) = 1577000.d0
    potion(333) = 1789600.d0;   potion(334) = 2005400.d0
    potion(335) = 2308000.d0;   potion(336) = 2536000.d0
    potion(337) = 2771000.d0;   potion(338) = 3250000.d0
    potion(339) = 3509900.d0;   potion(340) = 9144000.d0
    potion(341) = 9873000.d0;   potion(342) = 10649000.d0
    potion(343) = 11541000.d0;  potion(344) = 12398000.d0
    potion(345) = 13253000.d0;  potion(346) = 14427000.d0
    potion(347) = 15162200.d0;  potion(348) = 65659877.d0
    potion(349) = 69137430.d0
    ! Fe
    potion(351) =  63737.704d0; potion(352) = 130655.40d0
    potion(353) = 247220.d0;    potion(354) = 442900.d0
    potion(355) = 604900.d0;    potion(356) = 798370.d0
    potion(357) = 1008000.d0;   potion(358) = 1218380.d0
    potion(359) = 1884000.d0;   potion(360) = 2114000.d0
    potion(361) = 2346000.d0;   potion(362) = 2668000.d0
    potion(363) = 2912000.d0;   potion(364) = 3163000.d0
    potion(365) = 3680000.d0;   potion(366) = 3946570.d0
    potion(367) = 10184000.d0;  potion(368) = 10951000.d0
    potion(369) = 11770000.d0;  potion(370) = 12708000.d0
    potion(371) = 13607000.d0;  potion(372) = 14505000.d0
    potion(373) = 15731000.d0;  potion(374) = 16500160.d0
    potion(375) = 71204137.d0;  potion(376) = 74829550.d0
  end subroutine ionpots

end module mod_mklinelist
