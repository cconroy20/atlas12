MODULE synthe_module
! ============================================================================
!
!  MODULE synthe_module
!
!  PURPOSE
!  -------
!  F90 translation of the SYNTHE spectral synthesis package.
!  Contains all shared data and procedures for spectral synthesis.
!  Requires atlas12_modules.f90 (mod_atlas_data) for ATLAS library routines
!  (JOSH, READIN, KAPP, COMPUTE_ONE_POP, FREEFF) and atmosphere state.
!
!  SECTIONS (in order)
!  -------------------
!   1. Module-level shared data (replaces all COMMON blocks)
!   2. Utility routines: init_voigt_tables, voigt_profile, air_to_vac,
!      vac_to_air, interp_quadratic, trapz_integrate, parabolic_coeffs
!   3. Abundance/element data: load_abundances
!   4. Hydrogen oscillator strengths: hydrogen_oscillator_strength
!   5. Hydrogen line profiles: hydrogen_line_profile,
!      stark_quasistatic_profile, expint1
!   6. He I line profiles: he1_generic_profile, he1_4471_profile,
!      he1_4026_profile, he1_4387_profile, he1_4921_profile,
!      he1_griem_profile, he1_dimitri_profile, read_he1_stark_tables
!   7. Line opacity: compute_line_opacity
!   8. Atmosphere setup: run_xnfpelsyn (replaces standalone xnfpelsyn;
!      accesses mod_atlas_data via USE with renamed imports)
!
!  NOTES ON PRECISION
!  ------------------
!  The original code was deliberately mixed precision: REAL*4 for most
!  atmospheric quantities (for speed), REAL*8 for wavelengths and energy
!  levels.  That split is preserved here.  Variables that were explicitly
!  REAL*8 in the original remain REAL(8); everything else is REAL(4).
!
! ============================================================================

  USE mod_atlas_data, only: COMPUTE_ONE_POP, KAPP, FREEFF, DATADIR

  IMPLICIT NONE

  ! Public interface -- all module variables and procedures are public
  ! (removing the PRIVATE default since PROGRAM SYNTHE needs all shared data)

  ! ====================================================================
  !  MODULE-LEVEL SHARED DATA
  !  Replaces all COMMON blocks from the original F77 code.
  ! ====================================================================

  ! --- COMMON /ELEM/ ---
  ! Element abundance and atomic mass data, populated by load_abundances.
  !   abund(Z)  : log epsilon abundance relative to H = 0.911
  !   atmass(Z) : atomic mass (amu)
  !   elem(Z)   : 2-character element symbol
  REAL(4),          SAVE :: abund(99)
  REAL(4),          SAVE :: atmass(99)
  CHARACTER(LEN=2), SAVE :: elem(99)

  ! --- COMMON /H1TAB/ ---
  ! voigt_profile function look-up tables, populated by init_voigt_tables.
  !   H0TAB(i) = exp(-v^2)              at v = (i-1)/VSTEPS
  !   H1TAB(i) = H_1(v) first-order a correction kernel
  !   H2TAB(i) = (1 - 2v^2)*exp(-v^2)  second derivative term
  ! (Declared here; content moved out of separate module variables below)

  ! --- COMMON /EXTAB/ ---
  ! Fast exponential and E_1 tables, populated by init_voigt_tables.
  !   EXTAB(i)  = exp(-(i-1))
  !   EXTABF(i) = exp(-(i-1)*0.001)
  !   E1TAB(i)  = E_1(i*0.01)

  ! ====================================================================
  !  PHYSICAL CONSTANTS
  !  Centralised here so every subroutine uses the same literal values.
  !  All four forms of c are needed because different parts of the code
  !  work in different unit systems inherited from the F77 original.
  ! ====================================================================
  ! Constants also used by SPECTRV (keep in sync with standalone spectrv.f90)
  REAL(8), PARAMETER :: CLIGHT_NM_HZ     = 2.99792458D17  ! speed of light in nm·Hz
  REAL(8), PARAMETER :: PLANCK_PREFACTOR = 1.47439D-2      ! 2h/c² × (10¹⁵ Hz)³, CGS
  REAL(8), PARAMETER :: KMS_TO_CMS       = 1.0D5           ! km/s → cm/s

  REAL(8), PARAMETER :: CLIGHT_KMS = 299792.458D0    ! km/s
  REAL(8), PARAMETER :: CLIGHT_CMS = 2.99792458D10   ! cm/s
  REAL(8), PARAMETER :: CLIGHT_AHZ = 2.99792458D18   ! Å·Hz  (freq = CLIGHT_AHZ / lambda_Ang)

  ! ====================================================================
  !  LOOKUP TABLE SIZES
  !  Gathered here so the initialisation loops, init_voigt_tables call, and
  !  array declarations stay in sync if the resolution is ever changed.
  !    NTAB_EXP   : size of EXTAB / EXTABF   (integer + fractional parts)
  !    NTAB_E1    : size of E1TAB            (exponential integral E_1)
  !    NTAB_VOIGT : size of H0TAB/H1TAB/H2TAB  (voigt_profile profile kernels)
  ! ====================================================================
  INTEGER, PARAMETER :: NTAB_VOIGT = 2001

  ! ====================================================================
  !  voigt_profile APPROXIMATION THRESHOLD
  !  For damping parameter a < VOIGT_APPROX_LIMIT the first-order
  !  approximation H(0,a) ~ 1 - sqrt(pi)*a is used at line centre and
  !  H(v,a) ~ H0(v) + a*H1(v) for the wings, avoiding a full voigt_profile()
  !  call.  The threshold 0.2 gives ~1% accuracy (Hjerting 1938).
  !  The Boltzmann factor EXP(-ELO*HCKT) is computed via the intrinsic
  !  EXP() directly (replacing the former FASTEX table approximation).
  ! ====================================================================
  REAL(4), PARAMETER :: VOIGT_APPROX_LIMIT = 0.2

  ! --- COMMON /BHE/ ---
  ! He I bound-free and line opacity data for each depth point.
  INTEGER, PARAMETER :: kw = 99, mw = 139, mw6 = mw * 6
  REAL(4), SAVE :: BHE1(kw,29), AHE1(kw), SHE1(kw)
  REAL(4), SAVE :: BHE2(kw,6),  AHE2(kw), SHE2(kw)
  REAL(4), SAVE :: AHEMIN(kw),  SIGHE(kw)
  REAL(4), SAVE :: XNFPHE(kw,3), XNFHE(kw,2)

  ! --- COMMON /BHYD/ ---
  ! H I bound-free and line opacity data for each depth point.
  REAL(4), SAVE :: BHYD(kw,8), AHYD(kw), SHYD(kw), AH2P(kw)
  REAL(4), SAVE :: BMIN(kw), AHMIN(kw), SHMIN(kw)
  REAL(4), SAVE :: SIGH(kw), SIGH2(kw)
  REAL(4), SAVE :: AHLINE(kw), SHLINE(kw)
  REAL(4), SAVE :: XNFPH(kw,2), XNFH(kw)

  ! --- COMMON /H1TAB/ (declared as module variables) ---
  REAL(4), SAVE :: H0TAB(NTAB_VOIGT), H1TAB(NTAB_VOIGT), H2TAB(NTAB_VOIGT)

  ! --- COMMON /NLINES/ ---
  ! Wavelength grid parameters.
  REAL(8), SAVE :: WLBEG, WLEND, RESOLU, RATIO, RATIOLG, WBEGIN
  INTEGER, SAVE :: LENGTH, MLINES, IXWLBEG

  ! --- COMMON /RHOX/ ---
  REAL(4), SAVE :: RHOX(kw)
  INTEGER, SAVE :: NRHOX

  ! --- COMMON /STATE/ ---
  REAL(4), SAVE :: P(kw), XNE(kw), XNATOM(kw), RHO(kw), PTOTAL(kw)

  ! --- COMMON /TEMP/ ---
  REAL(4), SAVE :: T(kw), TKEV(kw), TK(kw), TLOG(kw)
  REAL(8), SAVE :: HKT(kw)   ! h*nu_0/(k*T) for stimulated-emission factor -- R8 for Boltzmann accuracy
  REAL(8), SAVE :: HCKT(kw)  ! h*c/(k*T) in cm -- exponent in EXP(-ELO*HCKT); R8 preserves high-excitation lines
  INTEGER, SAVE :: ITEMP

  ! --- COMMON /TURBPR/ ---
  REAL(4), SAVE :: VTURB(kw), PTURB(kw)
  REAL(4), SAVE :: TRBFDG, TRBCON, TRBPOW, TRBSND
  INTEGER, SAVE :: IFTURB

  ! --- COMMON /TXNXN/ ---
  REAL(8), SAVE :: EMERGE(kw)
  REAL(4), SAVE :: TXNXN(kw), BSTIM(kw), XNFH2(kw)

  ! --- COMMON /XNFDOP/ ---
  REAL(8), SAVE :: XNFPEL(mw6)  ! number density / rho for each ion (cm^{-3} g^{-1}) -- R8 for line strength precision
  REAL(8), SAVE :: DOPPLE(mw6)  ! Doppler width v/c for each ion -- R8 eliminates ~30 m/s round-off
  REAL(8), SAVE :: XNFDOP(mw6)  ! XNFPEL / DOPPLE: combined opacity factor -- R8 since derived from R8 inputs

  ! --- COMMON /LINDAT/ ---
  ! Full line data record.  REAL*8 for wavelengths/energy levels,
  ! REAL*4 for damping constants, INTEGER for quantum numbers.
  ! In PROGRAM SYNTHE the original code used
  !   EQUIVALENCE (LINDAT8(1),WL) and (LINDAT4(1),NELION)
  ! to read binary records; here we keep the named variables and the
  ! program fills them from LINDAT8/LINDAT4 after each READ.
  REAL(8),    SAVE :: WL, E, EP
  REAL(8),    SAVE :: LABEL(2), LABELP(2), OTHER1(2), OTHER2(2)
  REAL(8),    SAVE :: WLVAC, CENTER, CONCEN
  INTEGER,    SAVE :: NELION
  REAL(4),    SAVE :: GAMMAR, GAMMAS, GAMMAW, REF
  INTEGER,    SAVE :: NBLO, NBUP, ISO1, ISO2
  REAL(4),    SAVE :: X1, X2, GFLOG, XJ, XJP
  REAL(8),    SAVE :: CODE, ELO
  REAL(4),    SAVE :: GF, GS, GR, GW
  REAL(4),    SAVE :: DWL, DGFLOG, DGAMMAR, DGAMMAS, DGAMMAW
  REAL(4),    SAVE :: EXTRA1, EXTRA2, EXTRA3

  ! --- COMMON /BUFFER/ and /CONTIN/ ---
  ! Declared here so compute_line_opacity and SYNTHE share the same arrays.
  ! Sizes from the production PARAMETER statement.
  INTEGER, PARAMETER :: LENREC  = 1000000
  INTEGER, PARAMETER :: MAXLEN  = 100000001
  INTEGER, PARAMETER :: MAXPROF = 1000000
  INTEGER, PARAMETER :: MAXBUFF = MAXLEN + MAXPROF
  INTEGER, PARAMETER :: MAXLIN  = MAXBUFF + MAXPROF * 2
  REAL(4), SAVE :: BUFFER(MAXBUFF), PROFILE(MAXPROF)
  REAL(4), SAVE :: CONTINUUM(MAXBUFF)

  ! ====================================================================
  !  XNFPELSYN SHARED DATA
  !
  !  These arrays replace the fort.10 file written by XNFPELSYN and read
  !  by SYNTHE / SPECTRV.  They are filled by run_xnfpelsyn() and
  !  consumed by synthe_spectrv without any intermediate I/O.
  !
  !  Dimensions:
  !    kw   = 99   (max depth points)
  !    mw   = 139  (max elements + molecules)
  !    mm   = 100  (mw - 39; molecule count for IDMOL/MOMASS)
  !    nf   = 1131 (max continuum frequency points = 3*(nedge_max-1))
  !    me   = 377  (max continuum edge points)
  ! ====================================================================

  INTEGER, PARAMETER :: mm  = mw - 39   ! number of molecule entries
  INTEGER, PARAMETER :: nf  = 1131      ! max continuum frequency points
  INTEGER, PARAMETER :: me  = 377       ! max continuum edge points

  ! --- Edge frequency grid (written once; used by both SYNTHE and SPECTRV) ---
  REAL(8), SAVE :: frqedg_m(me)         ! frequency at each edge (Hz)
  REAL(8), SAVE :: wledge_m(me)         ! wavelength at each edge (nm)
  REAL(8), SAVE :: cmedge_m(me)         ! wavenumber at each edge (cm^{-1})
  REAL(8), SAVE :: idmol_m(mm)          ! molecule code identifiers
  REAL(8), SAVE :: momass_m(mm)         ! molecule masses
  INTEGER, SAVE :: nedge_m              ! number of edges

  ! --- Continuum frequency grid (for SPECTRV continuum interpolation) ---
  REAL(8), SAVE :: freqset_m(nf)        ! continuum frequency grid points
  INTEGER, SAVE :: numnu_m              ! number of grid points = 3*(nedge_m-1)

  ! --- Depth-structure arrays (replaces the depth-structure record on unit 10) ---
  ! These are filled by run_xnfpelsyn from the ATLAS COMMONs after READIN.
  ! Declared REAL(8) to match the original REAL*8 COMMON variables.
  REAL(8), SAVE :: xf_t(kw), xf_tkev(kw), xf_tk(kw), xf_hkt(kw)
  REAL(8), SAVE :: xf_tlog(kw), xf_hckt(kw)
  REAL(8), SAVE :: xf_p(kw), xf_xne(kw), xf_xnatom(kw), xf_rho(kw)
  REAL(8), SAVE :: xf_rhox(kw), xf_vturb(kw)
  REAL(8), SAVE :: xf_xnfh(kw), xf_xnfhe(kw,2), xf_xnfh2(kw)

  ! --- Per-depth continuum opacity tables ---
  ! continall_m(nu,j) = log10(kappa_abs + kappa_scat) at freq nu, depth j
  ! contabs_m(nu,j)   = log10(kappa_abs)
  ! contscat_m(nu,j)  = log10(kappa_scat)
  REAL(8), SAVE :: continall_m(nf, kw)
  REAL(8), SAVE :: contabs_m  (nf, kw)
  REAL(8), SAVE :: contscat_m (nf, kw)

  ! --- Per-depth ion populations and Doppler widths ---
  ! xnfpel_m(ion, elem, j) : number density of ion state (ion=1..6) of
  !                           element/molecule elem at depth j
  ! dopple_m(ion, elem, j) : Doppler width (fraction of c) for same
  REAL(8), SAVE :: xnfpel_m(6, mw, kw)
  REAL(8), SAVE :: dopple_m (6, mw, kw)

  ! ====================================================================
  !  SPECTRV SHARED INTERFACE DATA
  !
  !  These arrays replace the unit-9 / unit-8 I/O between SYNTHE and
  !  SPECTRV when the two programs are merged into synthe_spectrv.f90.
  !  They are declared here (in synthe_module) so that both the SYNTHE
  !  half and the SPECTRV half of the merged program can access them via
  !  USE synthe_module without passing large arrays through argument lists.
  !
  !  asynth_sv(j)          : net line opacity at depth j for the current
  !                          wavelength point (set by the SYNTHE transpose
  !                          section, consumed by the SPECTRV wavelength loop).
  !  lindat8_sv(14)        : REAL*8 line parameter record (wl, e, ep, labels,
  !                          wlvac, center, concen).
  !  lindat4_sv(28)        : REAL*4 line parameter record (nelion, damping, …).
  !  alinec_sv(kw)         : depth vector of line-centre opacities for the
  !                          current line identification record.
  !  nlines_sv             : number of line-centre identification records
  !                          (= n9 at the end of SYNTHE section 9).
  !
  !  SPECTRV continuum tables (filled once from unit 10 before the main loop):
  !  contabs_sv(3,medge,kw)  : log10 continuum absorption per (point,edge,depth)
  !  contscat_sv(3,medge,kw) : log10 continuum scattering
  !
  !  SPECTRV depth-dependent work arrays (filled once after READIN):
  !  bfudge_sv(kw)  : source-function fudge factor per depth
  !  fscat_sv(kw)   : scattering fraction per depth
  !  surf_sv(20)    : pure-continuum surface flux/intensity per angle quadrature
  !                   point (set by first JOSH call, used in residual computation)
  !
  !  SPECTRV scalar run parameters (filled from unit 25):
  !  rhoxj_sv, r1_sv, r101_sv, ph1_sv, pc1_sv, psi1_sv : see spectrv comments
  !  slope_sv, origin_sv : ASCII plot scale factors
  ! ====================================================================

  INTEGER, PARAMETER :: medge = 377   ! max continuum edge points (matches spectrv)

  REAL(4), SAVE :: asynth_sv(kw)
  REAL(8), SAVE :: lindat8_sv(14)
  REAL(4), SAVE :: lindat4_sv(28)
  REAL(4), SAVE :: alinec_sv(kw)
  INTEGER, SAVE :: nlines_sv

  REAL(8), SAVE :: contabs_sv(3, medge, kw)
  REAL(8), SAVE :: contscat_sv(3, medge, kw)

  REAL(8), SAVE :: bfudge_sv(kw)
  REAL(8), SAVE :: fscat_sv(kw)
  REAL(8), SAVE :: surf_sv(20)

  REAL(8), SAVE :: rhoxj_sv, r1_sv, r101_sv
  REAL(8), SAVE :: ph1_sv, pc1_sv, psi1_sv
  REAL(8), SAVE :: slope_sv, origin_sv

  ! Ground-state Boltzmann populations for bfudge computation.
  ! These replicate the old COMMON B-array ground-state values that
  ! the F77 PFSAHA used to store as a side effect.
  !   bhyd_gs(j)  = g_1 * exp(-E_1*HCKT(j)) / U(H I) * F(H I) * XNATOM * XABUND_H
  !   bc1_gs(j)   = same for C I ground state
  !   bc2_gs(j)   = same for C II ground state
  !   bsi1_gs(j)  = same for Si I ground state
  !   bsi2_gs(j)  = same for Si II ground state
  ! Filled by compute_bfudge_pops() in run_xnfpelsyn.
  REAL(8), SAVE :: bhyd_gs(kw)
  REAL(8), SAVE :: bc1_gs(kw), bc2_gs(kw)
  REAL(8), SAVE :: bsi1_gs(kw), bsi2_gs(kw)

CONTAINS

! ============================================================================
!  init_voigt_tables -- pre-tabulate the voigt_profile function on a uniform velocity grid.
!
!  Must be called once before voigt_profile is used.
!
!  Arguments:
!    VSTEPS  -- number of table steps per Doppler half-width (e.g. 200.)
!    N       -- number of table points (must equal NTAB_VOIGT = 2001)
!
!  The routine:
!    1. Populates H1TAB via quadratic interpolation (interp_quadratic) of the tabulated
!       H_1(v) kernel from Kurucz.
!    2. Computes H0TAB = exp(-v^2) and H2TAB = (1 - 2v^2)*exp(-v^2).
!
!  Note: EXTAB/EXTABF/E1TAB tables and FASTEX/FASTE1 have been removed;
!  the intrinsic EXP() is used directly throughout.
! ============================================================================
  SUBROUTINE init_voigt_tables(vsteps, n)
    REAL(4), INTENT(IN) :: vsteps
    INTEGER, INTENT(IN) :: n

    ! Tabulated H_1(v) kernel at irregular v spacing (from Kurucz SYNTHE)
    REAL(4), PARAMETER :: TABVI(81) = [ &
       0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  &
       1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  &
       2.0,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,  &
       3.0,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  &
       4.0,  4.2,  4.4,  4.6,  4.8,  5.0,  5.2,  5.4,  5.6,  5.8,  &
       6.0,  6.2,  6.4,  6.6,  6.8,  7.0,  7.2,  7.4,  7.6,  7.8,  &
       8.0,  8.2,  8.4,  8.6,  8.8,  9.0,  9.2,  9.4,  9.6,  9.8,  &
      10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8,  &
      12.0 ]

    REAL(4), PARAMETER :: TABH1(81) = [ &
      -1.12838,  -1.10596,  -1.04048,  -0.93703,  -0.80346,  -0.64945,  &
      -0.48552,  -0.32192,  -0.16772,  -0.03012,   0.08594,   0.17789,  &
       0.24537,   0.28981,   0.31394,   0.32130,   0.31573,   0.30094,  &
       0.28027,   0.25648,   0.231726,  0.207528,  0.184882,  0.164341, &
       0.146128,  0.130236,  0.116515,  0.104739,  0.094653,  0.086005, &
       0.078565,  0.072129,  0.066526,  0.061615,  0.057281,  0.053430, &
       0.049988,  0.046894,  0.044098,  0.041561,  0.039250,  0.035195, &
       0.031762,  0.028824,  0.026288,  0.024081,  0.022146,  0.020441, &
       0.018929,  0.017582,  0.016375,  0.015291,  0.014312,  0.013426, &
       0.012620,  0.0118860, 0.0112145, 0.0105990, 0.0100332, 0.0095119,&
       0.0090306, 0.0085852, 0.0081722, 0.0077885, 0.0074314, 0.0070985,&
       0.0067875, 0.0064967, 0.0062243, 0.0059688, 0.0057287, 0.0055030,&
       0.0052903, 0.0050898, 0.0049006, 0.0047217, 0.0045526, 0.0043924,&
       0.0042405, 0.0040964, 0.0039595 ]

    INTEGER :: i
    REAL(4) :: vv

    ! --- Build H0TAB as a v-grid scratch array for interp_quadratic, then overwrite ---
    DO i = 1, n
      H0TAB(i) = REAL(i - 1) / vsteps
    END DO

    ! Interpolate H_1 kernel onto the uniform v grid -> H1TAB
    CALL interp_quadratic(TABVI, TABH1, 81, H0TAB, H1TAB, n)

    ! Now fill H0TAB = exp(-v^2) and H2TAB = (1 - 2v^2)*exp(-v^2)
    DO i = 1, n
      vv       = (REAL(i - 1) / vsteps)**2
      H0TAB(i) = EXP(-vv)
      H2TAB(i) = H0TAB(i) * (1.0 - vv - vv)
    END DO

  END SUBROUTINE init_voigt_tables


! ============================================================================
!  voigt_profile -- fast approximation to the voigt_profile profile H(v, a).
!
!  H(v, a) is normalised so that integral H dv = sqrt(pi), i.e. the same
!  convention as used throughout SYNTHE / ATLAS.
!
!  Arguments:
!    V  -- dimensionless frequency offset in Doppler widths (>= 0)
!    A  -- damping parameter = gamma / (4*pi * Delta_nu_D)
!
!  Algorithm (Kurucz):
!    Three regimes:
!      a < VOIGT_APPROX_LIMIT, v <= 10:  quadratic Taylor expansion in a using pre-tabulated
!                          H0, H1, H2 (fast, accurate to ~a^3)
!      a < VOIGT_APPROX_LIMIT, v > 10:   far-wing approximation H ~ 0.5642 * a / v^2
!      a >= 0.2, a+v<=3.2: higher-order Taylor expansion with rational correction
!      a >= 0.2, large:    asymptotic Lorentz limit H ~ a * 0.79788 / (a^2+v^2)
!                           with correction terms
!
!  Requires: H0TAB, H1TAB, H2TAB initialised by init_voigt_tables.
! ============================================================================
  FUNCTION voigt_profile(v, a) RESULT(result)
    REAL(4), INTENT(IN) :: v, a
    REAL(4)             :: result
    INTEGER             :: iv
    REAL(4)             :: vv, aa, u, aau, vvu, uu
    REAL(4)             :: hh1, hh2, hh3, hh4

    ! Table index: 200 steps per Doppler width, offset 1.5 for rounding
    iv = INT(v * 200.0 + 1.5)

    IF (a >= VOIGT_APPROX_LIMIT) THEN

      IF (a > 1.4 .OR. a + v > 3.2) THEN
        ! --- Asymptotic / Lorentz regime ---
        aa     = a * a
        vv     = v * v
        u      = (aa + vv) * 1.4142
        result = a * 0.79788 / u
        IF (a > 100.0) RETURN
        aau    = aa / u
        vvu    = vv / u
        uu     = u * u
        result = (((((aau - 10.0*vvu) * aau * 3.0 + 15.0*vvu*vvu) &
                   + 3.0*vv - aa) / uu + 1.0)) * result
        RETURN
      END IF

      ! --- Intermediate regime: Taylor expansion in a with rational correction ---
      vv  = v * v
      hh1 = H1TAB(iv) + H0TAB(iv) * 1.12838
      hh2 = H2TAB(iv) + hh1 * 1.12838 - H0TAB(iv)
      hh3 = (1.0 - H2TAB(iv)) * 0.37613 - hh1 * 0.66667 * vv + hh2 * 1.12838
      hh4 = (3.0 * hh3 - hh1) * 0.37613 + H0TAB(iv) * 0.66667 * vv * vv
      result = ((((hh4*a + hh3)*a + hh2)*a + hh1)*a + H0TAB(iv)) * &
               (((-0.122727278*a + 0.532770573)*a - 0.96284325)*a + 0.979895032)
      RETURN

    END IF

    ! --- Small damping regime (a < VOIGT_APPROX_LIMIT) ---
    IF (v > 10.0) THEN
      ! Far wing: Lorentz / Gaussian transition
      result = 0.5642 * a / v**2
    ELSE
      ! Core: quadratic in a
      result = (H2TAB(iv) * a + H1TAB(iv)) * a + H0TAB(iv)
    END IF

  END FUNCTION voigt_profile


! ============================================================================
!  air_to_vac -- convert air wavelength to vacuum wavelength (both in nm).
!
!  Uses the Edlen dispersion formula for dry air at STP, iterated twice
!  for accuracy (the wavenumber argument is a vacuum wavenumber, so the
!  formula must be applied self-consistently).
!
!  Reference: Edlen (1953), with constants as used in SYNTHE/ATLAS.
! ============================================================================
  FUNCTION air_to_vac(w) RESULT(result)
    REAL(8), INTENT(IN) :: w       ! air wavelength (nm)
    REAL(8)             :: result
    REAL(8)             :: waven, wnew

    ! Edlen refractivity coefficients
    REAL(8), PARAMETER :: N0  = 1.0000834213D0
    REAL(8), PARAMETER :: NA  = 2406030.D0, CA = 1.30D10
    REAL(8), PARAMETER :: NB  =   15997.D0, CB = 3.89D9

    ! First iteration (use air wavenumber as approximation)
    waven = 1.D7 / w
    wnew  = w * (N0 + NA/(CA - waven**2) + NB/(CB - waven**2))

    ! Second iteration
    waven = 1.D7 / wnew
    wnew  = w * (N0 + NA/(CA - waven**2) + NB/(CB - waven**2))

    ! Third iteration (final)
    waven  = 1.D7 / wnew
    result = w * (N0 + NA/(CA - waven**2) + NB/(CB - waven**2))

  END FUNCTION air_to_vac


! ============================================================================
!  vac_to_air -- convert vacuum wavelength to air wavelength (both in nm).
!
!  Single application of the Edlen formula (no iteration needed since
!  the vacuum wavenumber is exact).
! ============================================================================
  FUNCTION vac_to_air(w) RESULT(result)
    REAL(8), INTENT(IN) :: w       ! vacuum wavelength (nm)
    REAL(8)             :: result
    REAL(8)             :: waven

    REAL(8), PARAMETER :: N0  = 1.0000834213D0
    REAL(8), PARAMETER :: NA  = 2406030.D0, CA = 1.30D10
    REAL(8), PARAMETER :: NB  =   15997.D0, CB = 3.89D9

    waven  = 1.D7 / w
    result = w / (N0 + NA/(CA - waven**2) + NB/(CB - waven**2))

  END FUNCTION vac_to_air


! ============================================================================
!  interp_quadratic -- remap a function from an old grid to a new grid using piecewise
!          quadratic (parabolic) interpolation with blended coefficients.
!
!  Arguments:
!    XOLD(NOLD)  -- old x grid (must be monotonically increasing)
!    FOLD(NOLD)  -- function values on old grid
!    NOLD        -- number of old grid points
!    XNEW(NNEW)  -- new x grid points (must be monotonically increasing)
!    FNEW(NNEW)  -- output: interpolated function values on new grid
!    NNEW        -- number of new grid points
!
!  The algorithm fits piecewise quadratics f = A + B*x + C*x^2 to each
!  interval, blending the "forward" and "backward" parabola fits at each
!  point to reduce Runge-like oscillations.
!
!  This routine is used by init_voigt_tables to map the irregular H_1(v) data
!  onto the uniform velocity grid.
! ============================================================================
  SUBROUTINE interp_quadratic(xold, fold, nold, xnew, fnew, nnew)
    REAL(4), INTENT(IN)  :: xold(nold), fold(nold)
    INTEGER, INTENT(IN)  :: nold, nnew
    REAL(4), INTENT(IN)  :: xnew(nnew)
    REAL(4), INTENT(OUT) :: fnew(nnew)

    INTEGER :: k, l, l1, l2, ll
    REAL(4) :: a, b, c, d, wt
    REAL(4) :: abac, bbac, cbac, afor, bfor, cfor

    l  = 2
    ll = 0
    a  = 0.0;  b  = 0.0;  c  = 0.0
    abac = 0.0; bbac = 0.0; cbac = 0.0
    afor = 0.0; bfor = 0.0; cfor = 0.0

    DO k = 1, nnew

      ! Advance l until xnew(k) < xold(l)
      DO WHILE (xnew(k) >= xold(l) .AND. l <= nold)
        l = l + 1
        IF (l > nold) EXIT
      END DO

      IF (l > nold) THEN
        ! Extrapolate beyond right end: use last linear segment
        IF (l /= ll) THEN
          l  = MIN(nold, l)
          c  = 0.0
          b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
          a  = fold(l) - xold(l) * b
          ll = l
        END IF

      ELSE IF (l == ll) THEN
        ! Same interval as last point: reuse coefficients
        CONTINUE

      ELSE IF (l == 2) THEN
        ! Left end: use linear segment
        c  = 0.0
        b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
        a  = fold(l) - xold(l) * b
        ll = l

      ELSE
        l1 = l - 1
        l2 = l - 2

        ! Compute backward parabola (through points l2, l1, l)
        IF (l > ll + 1 .OR. l == 3) THEN
          d    = (fold(l1) - fold(l2)) / (xold(l1) - xold(l2))
          cbac = fold(l)  / ((xold(l) - xold(l1)) * (xold(l)  - xold(l2))) + &
                (fold(l2) / (xold(l)  - xold(l2)) - &
                 fold(l1) / (xold(l)  - xold(l1))) / (xold(l1) - xold(l2))
          bbac = d - (xold(l1) + xold(l2)) * cbac
          abac = fold(l2) - xold(l2)*d + xold(l1)*xold(l2)*cbac
        ELSE
          ! Reuse forward as backward (already advanced by one)
          cbac = cfor
          bbac = bfor
          abac = afor
        END IF

        IF (l >= nold) THEN
          ! At or past right end: use backward parabola
          a  = abac;  b = bbac;  c = cbac
          ll = l
        ELSE
          ! Compute forward parabola (through points l1, l, l+1)
          d    = (fold(l) - fold(l1)) / (xold(l) - xold(l1))
          cfor = fold(l+1) / ((xold(l+1) - xold(l)) * (xold(l+1) - xold(l1))) + &
                (fold(l1) / (xold(l+1) - xold(l1)) - &
                 fold(l)  / (xold(l+1) - xold(l)))  / (xold(l) - xold(l1))
          bfor = d - (xold(l) + xold(l1)) * cfor
          afor = fold(l1) - xold(l1)*d + xold(l)*xold(l1)*cfor

          ! Blend forward and backward based on relative curvature
          wt = 0.0
          IF (ABS(cfor) > 0.0) wt = ABS(cfor) / (ABS(cfor) + ABS(cbac))
          a  = afor + wt * (abac - afor)
          b  = bfor + wt * (bbac - bfor)
          c  = cfor + wt * (cbac - cfor)
          ll = l
        END IF

      END IF

      fnew(k) = a + (b + c * xnew(k)) * xnew(k)

    END DO

  END SUBROUTINE interp_quadratic


! ============================================================================
!  trapz_integrate -- compute the definite integral of a piecewise-quadratic fit to f(x).
!
!  Arguments:
!    X(N)    -- independent variable grid (monotonically increasing)
!    F(N)    -- function values at grid points
!    FINT    -- output: integral from X(1) to X(N), plus START
!    N       -- number of points
!    START   -- initial value added to the integral (e.g. 0.)
!
!  Algorithm: fits piecewise quadratics via parabolic_coeffs, then analytically
!  integrates each segment.
! ============================================================================
  SUBROUTINE trapz_integrate(x, f, fint, n, start)
    INTEGER, INTENT(IN)  :: n
    REAL(4), INTENT(IN)  :: x(n), f(n), start
    REAL(4), INTENT(OUT) :: fint

    REAL(4) :: a(1000), b(1000), c(1000)
    INTEGER :: i

    CALL parabolic_coeffs(f, x, a, b, c, n)
    fint = start

    DO i = 1, n - 1
      fint = fint + (a(i) + b(i)/2.0*(x(i+1) + x(i)) + &
             c(i)/3.0*((x(i+1) + x(i))*x(i+1) + x(i)*x(i))) * (x(i+1) - x(i))
    END DO

  END SUBROUTINE trapz_integrate


! ============================================================================
!  parabolic_coeffs -- compute piecewise-quadratic polynomial coefficients for trapz_integrate.
!
!  For each interval [X(i), X(i+1)], fits f(x) ~ A(i) + B(i)*x + C(i)*x^2
!  using the neighbouring points, then blends adjacent fits to reduce
!  Runge oscillations.
!
!  Arguments:
!    F(N), X(N)  -- input function and grid
!    A(N), B(N), C(N)  -- output polynomial coefficients
!    N           -- number of points
! ============================================================================
  SUBROUTINE parabolic_coeffs(f, x, a, b, c, n)
    INTEGER, INTENT(IN)  :: n
    REAL(4), INTENT(IN)  :: f(n), x(n)
    REAL(4), INTENT(OUT) :: a(n), b(n), c(n)

    INTEGER :: j, j1
    REAL(4) :: d, wt
    INTEGER :: n1

    ! Endpoints: linear fit (zero curvature)
    c(1) = 0.0
    b(1) = (f(2) - f(1)) / (x(2) - x(1))
    a(1) = f(1) - x(1) * b(1)

    n1   = n - 1
    c(n) = 0.0
    b(n) = (f(n) - f(n1)) / (x(n) - x(n1))
    a(n) = f(n) - x(n) * b(n)

    IF (n == 2) RETURN

    ! Interior points: quadratic fit through each triplet
    DO j = 2, n1
      j1   = j - 1
      d    = (f(j) - f(j1)) / (x(j) - x(j1))
      c(j) = f(j+1) / ((x(j+1) - x(j)) * (x(j+1) - x(j1))) - &
             f(j)   / ((x(j)   - x(j1)) * (x(j+1) - x(j))) + &
             f(j1)  / ((x(j)   - x(j1)) * (x(j+1) - x(j1)))
      b(j) = d - (x(j) + x(j1)) * c(j)
      a(j) = f(j1) - x(j1)*d + x(j)*x(j1)*c(j)
    END DO

    ! Force linear fit at points 2 and 3 (suppress edge oscillations)
    c(2) = 0.0
    b(2) = (f(3) - f(2)) / (x(3) - x(2))
    a(2) = f(2) - x(2) * b(2)
    c(3) = 0.0
    b(3) = (f(4) - f(3)) / (x(4) - x(3))
    a(3) = f(3) - x(3) * b(3)

    ! Blend adjacent quadratics to suppress oscillations
    DO j = 2, n1
      IF (ABS(c(j)) < TINY(c(j))) CYCLE
      j1 = j + 1
      wt   = ABS(c(j1)) / (ABS(c(j1)) + ABS(c(j)))
      a(j) = a(j1) + wt * (a(j) - a(j1))
      b(j) = b(j1) + wt * (b(j) - b(j1))
      c(j) = c(j1) + wt * (c(j) - c(j1))
    END DO

    a(n1) = a(n)
    b(n1) = b(n)
    c(n1) = c(n)

  END SUBROUTINE parabolic_coeffs

! ============================================================================
!  load_abundances -- initialise element abundances, atomic masses, and symbols.
!
!  Replaces the original SUBROUTINE load_abundances which used COMMON /ELEM/.
!  Data source: Anders & Grevesse (1989), Geochimica et Cosmochimica Acta,
!  53, 197-214, as compiled by Kurucz.
!
!  Abundances are log epsilon relative to a hydrogen abundance of 0.911
!  (i.e. H is set to -0.04 rather than the conventional 12.00).
!  Elements with no known/relevant abundance are set to -20.00.
!
!  This subroutine must be called once before the abundance array is used.
! ============================================================================
  SUBROUTINE load_abundances()

    !                      1H       2He
    abund(1:2) = [ 0.911,  0.089 ]

    !                3Li      4Be      5B       6C       7N       8O
    abund(3:10) = [ -10.88, -10.89,  -9.44,   -3.48,   -3.99,   -3.11, &
    !                9F      10Ne
                     -7.48,  -3.95 ]

    !               11Na     12Mg     13Al     14Si     15P      16S
    abund(11:18) = [ -5.71,  -4.46,  -5.57,   -4.49,   -6.59,   -4.83, &
    !               17Cl     18Ar
                     -6.54,  -5.48 ]

    !               19K      20Ca     21Sc     22Ti     23V      24Cr
    abund(19:26) = [ -6.92,  -5.68,  -8.94,   -7.05,   -8.04,   -6.37, &
    !               25Mn     26Fe
                     -6.65,  -4.37 ]

    !               27Co     28Ni     29Cu     30Zn     31Ga     32Ge
    abund(27:34) = [ -7.12,  -5.79,  -7.83,   -7.44,   -9.16,   -8.63, &
    !               33As     34Se
                     -9.67,  -8.69 ]

    !               35Br     36Kr     37Rb     38Sr     39Y      40Zr
    abund(35:42) = [ -9.41,  -8.81,  -9.44,   -9.14,   -9.80,   -9.44, &
    !               41Nb     42Mo
                    -10.62, -10.12 ]

    !               43Tc     44Ru     45Rh     46Pd     47Ag     48Cd
    abund(43:50) = [-20.00, -10.20, -10.92,  -10.35,  -11.10,  -10.18, &
    !               49In     50Sn
                    -10.38, -10.04 ]

    !               51Sb     52Te     53I      54Xe     55Cs     56Ba
    abund(51:58) = [-11.04,  -9.80, -10.53,   -9.81,  -10.92,   -9.91, &
    !               57La     58Ce
                    -10.82, -10.49 ]

    !               59Pr     60Nd     61Pm     62Sm     63Eu     64Gd
    abund(59:66) = [-11.33, -10.54, -20.00,  -11.04,  -11.53,  -10.92, &
    !               65Tb     66Dy
                    -12.14, -10.94 ]

    !               67Ho     68Er     69Tm     70Yb     71Lu     72Hf
    abund(67:74) = [-11.78, -11.11, -12.04,  -10.96,  -11.28,  -11.16, &
    !               73Ta     74W
                    -11.91, -10.93 ]

    !               75Re     76Os     77Ir     78Pt     79Au     80Hg
    abund(75:82) = [-11.77, -10.59, -10.69,  -10.24,  -11.03,  -10.95, &
    !               81Tl     82Pb
                    -11.14, -10.19 ]

    !               83Bi     84Po     85At     86Rn     87Fr     88Ra
    abund(83:90) = [-11.33, -20.00, -20.00,  -20.00,  -20.00,  -20.00, &
    !               89Ac     90Th
                    -20.00, -11.92 ]

    !               91Pa     92U      93Np     94Pu     95Am     96Cm
    abund(91:99) = [-20.00, -12.51, -20.00,  -20.00,  -20.00,  -20.00, &
    !               97Bk     98Cf     99Es
                    -20.00, -20.00, -20.00 ]

    ! -------------------------------------------------------------------
    !  Atomic masses (amu) -- Anders & Grevesse (1989) / Kurucz
    ! -------------------------------------------------------------------
    !            1H      2He
    atmass(1:2) = [ 1.008,  4.003 ]

    atmass(3:12)  = [  6.939,  9.013, 10.81,  12.01,  14.01,  16.00, &
                       19.00, 20.18,  22.99,  24.31 ]
    atmass(13:22) = [ 26.98,  28.09,  30.98,  32.07,  35.45,  39.95, &
                      39.10,  40.08,  44.96,  47.90 ]
    atmass(23:32) = [ 50.94,  52.00,  54.94,  55.85,  58.94,  58.71, &
                      63.55,  65.37,  69.72,  72.60 ]
    atmass(33:42) = [ 74.92,  78.96,  79.91,  83.80,  85.48,  87.63, &
                      88.91,  91.22,  92.91,  95.95 ]
    atmass(43:52) = [ 99.00, 101.1,  102.9,  106.4,  107.9,  112.4, &
                     114.8,  118.7,  121.8,  127.6 ]
    atmass(53:62) = [126.9,  131.3,  132.9,  137.4,  138.9,  140.1, &
                     140.9,  144.3,  147.0,  150.4 ]
    atmass(63:72) = [152.0,  157.3,  158.9,  162.5,  164.9,  167.3, &
                     168.9,  173.0,  175.0,  178.5 ]
    atmass(73:82) = [181.0,  183.9,  186.3,  190.2,  192.2,  195.1, &
                     197.0,  200.6,  204.4,  207.2 ]
    atmass(83:92) = [209.0,  210.0,  211.0,  222.0,  223.0,  226.1, &
                     227.1,  232.0,  231.0,  238.0 ]
    atmass(93:99) = [237.0,  244.0,  243.0,  247.0,  247.0,  251.0,  254.0 ]

    ! -------------------------------------------------------------------
    !  Element symbols (2-character, left-justified)
    !  Replaces the original Hollerith DATA ELEM/ 2HH , 2HHE, ... /
    !  stored in a REAL*8 array -- a notably fragile F77 idiom.
    ! -------------------------------------------------------------------
    elem(1:10)  = [ 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE' ]
    elem(11:20) = [ 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA' ]
    elem(21:30) = [ 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN' ]
    elem(31:40) = [ 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR' ]
    elem(41:50) = [ 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN' ]
    elem(51:60) = [ 'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND' ]
    elem(61:70) = [ 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB' ]
    elem(71:80) = [ 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG' ]
    elem(81:90) = [ 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH' ]
    elem(91:99) = [ 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES' ]

  END SUBROUTINE load_abundances


! ============================================================================
!  hydrogen_oscillator_strength -- hydrogen oscillator strength f_{nm} for transition n -> m.
! ============================================================================
  FUNCTION hydrogen_oscillator_strength(n, m) RESULT(result)
    INTEGER, INTENT(IN) :: n, m
    REAL(4)             :: result

    INTEGER, SAVE :: nstr = 0, mstr = 0
    REAL(4),  SAVE :: fnm = 0.0

    REAL(4) :: xn, xm, xmn, xmn12, fkn, ginf, gca, fk, wtc, wt

    result = 0.0
    IF (m <= n) RETURN

    IF (n /= nstr) THEN
      xn   = REAL(n)
      ginf = 0.2027 / xn**0.71
      gca  = 0.124  / xn
      fkn  = xn * 1.9603
      wtc  = 0.45 - 2.4 / xn**3 * (xn - 1.0)
      nstr = n
    END IF

    IF (m /= mstr) THEN
      xm    = REAL(m)
      xmn   = REAL(m - n)
      fk    = fkn * (xm / (xmn * (xm + REAL(n))))**3
      xmn12 = xmn**1.2
      wt    = (xmn12 - 1.0) / (xmn12 + wtc)
      fnm   = fk * (1.0 - wt*ginf - (0.222 + gca/xm)*(1.0 - wt))
      mstr  = m
    END IF

    result = fnm

  END FUNCTION hydrogen_oscillator_strength


! ============================================================================
!  expint1 -- evaluation of the exponential integral E_1(x) for x >= 0.
!
!  Used inside hydrogen_line_profile for the Holtsmark-Stark broadening calculation.
!
!  Two-branch implementation:
!    0 < x <= 1 : Abramowitz & Stegun 5.1.56 five-term polynomial in x,
!                 max relative error ~1.5e-7.  Handles x -> 0 correctly;
!                 the x<=0.01 one-term branch formerly here was less accurate
!                 than this polynomial and has been removed.
!    1 < x <= 30: Pade approximant for e^x * x * E1(x), max error ~5e-5.
!    x > 30     : returns 0 (E1(30) ~ 3e-15, negligible for Stark profiles).
!
!  The 5e-5 Pade error is ~200x smaller than the ~1% accuracy of the Stark
!  profile tables (stark_quasistatic_profile), so no further improvement is warranted.
! ============================================================================
  FUNCTION expint1(x) RESULT(result)
    REAL(4), INTENT(IN) :: x
    REAL(4)             :: result

    result = 0.0
    IF (x <= 0.0) RETURN

    IF (x <= 1.0) THEN
      ! A&S 5.1.56: E_1(x) = -ln(x) + a0 + a1*x + a2*x^2 + ... + a5*x^5
      ! Accurate to 2e-7 absolute; equivalent to the series -ln(x) - gamma + sum_n
      result = -LOG(x) - 0.57721566 + x * (0.99999193 + x * (-0.24991055 + x * &
               (0.05519968 + x * (-0.00976004 + x * 0.00107857))))
      RETURN
    END IF

    IF (x <= 30.0) THEN
      ! Pade approximant for e^x * x * E_1(x), from Abramowitz & Stegun 5.1.56
      result = (x*(x + 2.334733) + 0.25062) / (x*(x + 3.330657) + 1.681534) &
               / x * EXP(-x)
    END IF
    ! x > 30: E_1(x) ~ 3e-15 or smaller -- return 0

  END FUNCTION expint1


! ============================================================================
!  stark_quasistatic_profile -- quasistatic Stark profile S(beta, p) for hydrogen lines.
!
!  Generates the quasistatic ion field profile S(beta, p) used in hydrogen_line_profile.
!  The alpha and beta lines of the first three Lyman/Balmer/Paschen series
!  are treated explicitly; H18 profile is used for the rest.
!
!  Arguments:
!    B  -- reduced field strength beta = |Delta_nu| / F_0 * D_beta
!    P  -- plasma parameter (ratio of interparticle distance to Debye length)
!    N  -- lower principal quantum number
!    M  -- upper principal quantum number
!
!  Data source: tabulated corrections PROPBM from Kurucz (1993),
!  fitted to Vidal, Cooper & Smith (1973) profiles.
! ============================================================================
  FUNCTION stark_quasistatic_profile(b, p, n, m) RESULT(result)
    REAL(4), INTENT(IN) :: b, p
    INTEGER, INTENT(IN) :: n, m
    REAL(4)             :: result

    ! Interpolation grids
    REAL(4), PARAMETER :: PP(5)   = [ 0.0, 0.2, 0.4, 0.6, 0.8 ]
    REAL(4), PARAMETER :: BETA(15) = [ 1.0, 1.259, 1.585, 1.995, 2.512, 3.162, &
                                        3.981, 5.012, 6.310, 7.943, 10.0, 12.59, &
                                        15.85, 19.95, 25.12 ]

    ! Correction table PROPBM(IP, IBETA, ILINE), stored as 7 flat arrays of 75
    REAL(4), PARAMETER :: PROB1(75) = [ &
      -.980,-.967,-.948,-.918,-.873,-.968,-.949,-.921,-.879,-.821, &
      -.950,-.922,-.883,-.830,-.764,-.922,-.881,-.830,-.770,-.706, &
      -.877,-.823,-.763,-.706,-.660,-.806,-.741,-.682,-.640,-.625, &
      -.691,-.628,-.588,-.577,-.599,-.511,-.482,-.484,-.514,-.568, &
      -.265,-.318,-.382,-.455,-.531,-.013,-.167,-.292,-.394,-.478, &
       .166,-.056,-.216,-.332,-.415, .251, .035,-.122,-.237,-.320, &
       .221, .059,-.068,-.168,-.247, .160, .055,-.037,-.118,-.189, &
       .110, .043,-.022,-.085,-.147 ]  ! Lyman alpha

    REAL(4), PARAMETER :: PROB2(75) = [ &
      -.242, .060, .379, .671, .894,  .022, .314, .569, .746, .818, &
       .273, .473, .605, .651, .607,  .432, .484, .489, .442, .343, &
       .434, .366, .294, .204, .091,  .304, .184, .079,-.025,-.135, &
       .167, .035,-.082,-.189,-.290,  .085,-.061,-.183,-.287,-.374, &
       .032,-.127,-.249,-.344,-.418, -.024,-.167,-.275,-.357,-.420, &
      -.061,-.170,-.257,-.327,-.384, -.047,-.124,-.192,-.252,-.306, &
      -.043,-.092,-.142,-.190,-.238, -.038,-.070,-.107,-.146,-.187, &
      -.030,-.049,-.075,-.106,-.140 ]  ! Lyman beta

    REAL(4), PARAMETER :: PROB3(75) = [ &
      -.484,-.336,-.206,-.111,-.058, -.364,-.264,-.192,-.154,-.144, &
      -.299,-.268,-.250,-.244,-.246, -.319,-.333,-.337,-.336,-.337, &
      -.397,-.414,-.415,-.413,-.420, -.456,-.455,-.451,-.456,-.478, &
      -.446,-.441,-.446,-.469,-.512, -.358,-.381,-.415,-.463,-.522, &
      -.214,-.288,-.360,-.432,-.503, -.063,-.196,-.304,-.394,-.468, &
       .063,-.108,-.237,-.334,-.409,  .151,-.019,-.148,-.245,-.319, &
       .149, .016,-.091,-.177,-.246,  .115, .023,-.056,-.126,-.189, &
       .078, .021,-.036,-.091,-.145 ]  ! Balmer alpha

    REAL(4), PARAMETER :: PROB4(75) = [ &
      -.082, .163, .417, .649, .829,  .096, .316, .515, .660, .729, &
       .242, .393, .505, .556, .534,  .320, .373, .394, .369, .290, &
       .308, .274, .226, .152, .048,  .232, .141, .052,-.046,-.154, &
       .148, .020,-.094,-.200,-.299,  .083,-.070,-.195,-.299,-.385, &
       .031,-.130,-.253,-.348,-.422, -.023,-.167,-.276,-.359,-.423, &
      -.053,-.165,-.254,-.326,-.384, -.038,-.119,-.190,-.251,-.306, &
      -.034,-.088,-.140,-.190,-.239, -.032,-.066,-.103,-.144,-.186, &
      -.027,-.048,-.075,-.106,-.142 ]  ! Balmer beta

    REAL(4), PARAMETER :: PROB5(75) = [ &
      -.819,-.759,-.689,-.612,-.529, -.770,-.707,-.638,-.567,-.498, &
      -.721,-.659,-.595,-.537,-.488, -.671,-.617,-.566,-.524,-.497, &
      -.622,-.582,-.547,-.523,-.516, -.570,-.545,-.526,-.521,-.537, &
      -.503,-.495,-.496,-.514,-.551, -.397,-.418,-.448,-.492,-.547, &
      -.246,-.315,-.384,-.453,-.522, -.080,-.210,-.316,-.406,-.481, &
       .068,-.107,-.239,-.340,-.418,  .177,-.006,-.143,-.246,-.324, &
       .184, .035,-.082,-.174,-.249,  .146, .042,-.046,-.123,-.190, &
       .103, .036,-.027,-.088,-.146 ]  ! Paschen alpha

    REAL(4), PARAMETER :: PROB6(75) = [ &
      -.073, .169, .415, .636, .809,  .102, .311, .499, .639, .710, &
       .232, .372, .479, .531, .514,  .294, .349, .374, .354, .279, &
       .278, .253, .212, .142, .040,  .215, .130, .044,-.051,-.158, &
       .141, .015,-.097,-.202,-.300,  .080,-.072,-.196,-.299,-.385, &
       .029,-.130,-.252,-.347,-.421, -.022,-.166,-.275,-.359,-.423, &
      -.050,-.164,-.253,-.325,-.384, -.035,-.118,-.189,-.252,-.306, &
      -.032,-.087,-.139,-.190,-.240, -.029,-.064,-.102,-.143,-.185, &
      -.025,-.046,-.074,-.106,-.142 ]  ! Paschen beta

    REAL(4), PARAMETER :: PROB7(75) = [ &
       .005, .128, .260, .389, .504,  .004, .109, .220, .318, .389, &
      -.007, .079, .162, .222, .244, -.018, .041, .089, .106, .080, &
      -.026,-.003, .003,-.023,-.086, -.025,-.048,-.087,-.148,-.234, &
      -.008,-.085,-.165,-.251,-.343,  .018,-.111,-.223,-.321,-.407, &
       .032,-.130,-.255,-.354,-.431,  .014,-.148,-.269,-.359,-.427, &
      -.005,-.140,-.243,-.323,-.386,  .005,-.095,-.178,-.248,-.307, &
      -.002,-.068,-.129,-.187,-.241, -.007,-.049,-.094,-.139,-.186, &
      -.010,-.036,-.067,-.103,-.143 ]  ! Balmer 18 (generic high-n)

    ! Asymptotic correction coefficients C(IP, ILINE) and D(IP, ILINE)
    REAL(4), PARAMETER :: CC(5,7) = RESHAPE( [ &
      -18.396, 84.674, -96.273,   3.927,  55.191, &   ! Ly alpha
       95.740, 18.489,  14.902,  24.466,  42.456, &   ! Ly beta
      -25.088,145.882, -50.165,   7.902,  51.003, &   ! Ba alpha
       93.783, 10.066,   9.224,  20.685,  40.136, &   ! Ba beta
      -19.819, 94.981, -79.606,   3.159,  52.106, &   ! Pa alpha
      111.107, 11.910,   9.857,  21.371,  41.006, &   ! Pa beta
      511.318,  1.532,   4.044,  19.266,  41.812  &   ! Ba 18
      ], [5, 7] )

    REAL(4), PARAMETER :: DD(5,7) = RESHAPE( [ &
       11.801,  9.079,  -0.651, -11.071, -26.545, &
       -6.665, -7.136, -10.605,-15.882, -23.632, &
        7.872,  5.592,  -2.716, -12.180, -25.661, &
       -5.918, -6.501, -10.130,-15.588, -23.570, &
       10.938,  8.028,  -1.267, -11.375, -26.047, &
       -5.899, -6.381, -10.044,-15.574, -23.644, &
       -6.070, -4.528,  -8.759, -14.984, -23.956  &
      ], [5, 7] )

    INTEGER :: indx, im, ip, jm, jp, j, mmn
    REAL(4) :: b2, sb, wtpp, wtpm, wtbp, wtbm
    REAL(4) :: cbp, cbm, corr, pr1, pr2, wt, ent
    REAL(4) :: prop_im_jm, prop_ip_jm, prop_im_jp, prop_ip_jp
    REAL(4) :: c_val, d_val

    b2 = b * b
    sb = SQRT(b)

    IF (b > 500.0) THEN
      ! Pure asymptotic
      result = (1.5/sb + 27.0/b2) / b2
      RETURN
    END IF

    ! Select which line profile table to use
    mmn  = m - n
    indx = 7   ! default: Balmer-18
    IF (n <= 3 .AND. mmn <= 2) indx = 2*(n-1) + mmn

    ! Plasma parameter interpolation weights
    im   = MIN(INT(5.0*p) + 1, 4)
    ip   = im + 1
    wtpp = 5.0 * (p - PP(im))
    wtpm = 1.0 - wtpp

    IF (b > 25.12) THEN
      ! Asymptotic regime with correction
      c_val  = CC(ip,indx)*wtpp + CC(im,indx)*wtpm
      d_val  = DD(ip,indx)*wtpp + DD(im,indx)*wtpm
      corr   = 1.0 + d_val / (c_val + b*sb)
      result = (1.5/sb + 27.0/b2) / b2 * corr
      RETURN
    END IF

    ! Find beta bracket
    jp = 2
    DO j = 2, 15
      jp = j
      IF (b <= BETA(j)) EXIT
    END DO
    jm   = jp - 1
    wtbp = (b - BETA(jm)) / (BETA(jp) - BETA(jm))
    wtbm = 1.0 - wtbp

    ! Look up correction table for the four surrounding points
    ! PROPBM is logically PROPBM(ip, jp, indx); stored as PROB arrays of 75
    ! index = (ibeta-1)*5 + ip  within each PROB block
    CALL stark_profile_interp(indx, im, ip, jm, jp, &
                       PROB1, PROB2, PROB3, PROB4, PROB5, PROB6, PROB7, &
                       prop_im_jm, prop_ip_jm, prop_im_jp, prop_ip_jp)

    cbp  = prop_ip_jp*wtpp + prop_im_jp*wtpm
    cbm  = prop_ip_jm*wtpp + prop_im_jm*wtpm
    corr = 1.0 + cbp*wtbp + cbm*wtbm

    ! Inner approximate profile (transition region)
    pr1 = 0.0;  pr2 = 0.0
    wt  = MAX(MIN(0.5*(10.0 - b), 1.0), 0.0)
    IF (b <= 10.0) pr1 = 8.0 / (83.0 + (2.0 + 0.95*b2)*b)
    IF (b >= 8.0)  pr2 = (1.5/sb + 27.0/b2) / b2
    result = (pr1*wt + pr2*(1.0 - wt)) * corr

  END FUNCTION stark_quasistatic_profile

  ! Helper: extract four PROPBM values from the flat PROB arrays
  PURE SUBROUTINE stark_profile_interp(indx, im, ip, jm, jp, &
                                 P1, P2, P3, P4, P5, P6, P7, &
                                 vim_jm, vip_jm, vim_jp, vip_jp)
    INTEGER, INTENT(IN) :: indx, im, ip, jm, jp
    REAL(4), INTENT(IN) :: P1(75),P2(75),P3(75),P4(75),P5(75),P6(75),P7(75)
    REAL(4), INTENT(OUT):: vim_jm, vip_jm, vim_jp, vip_jp
    ! PROPBM(ip, ibeta, indx): index within a 75-element block = (ibeta-1)*5 + ip
    INTEGER :: i_im_jm, i_ip_jm, i_im_jp, i_ip_jp
    i_im_jm = (jm-1)*5 + im
    i_ip_jm = (jm-1)*5 + ip
    i_im_jp = (jp-1)*5 + im
    i_ip_jp = (jp-1)*5 + ip
    SELECT CASE (indx)
      CASE (1); vim_jm=P1(i_im_jm); vip_jm=P1(i_ip_jm); vim_jp=P1(i_im_jp); vip_jp=P1(i_ip_jp)
      CASE (2); vim_jm=P2(i_im_jm); vip_jm=P2(i_ip_jm); vim_jp=P2(i_im_jp); vip_jp=P2(i_ip_jp)
      CASE (3); vim_jm=P3(i_im_jm); vip_jm=P3(i_ip_jm); vim_jp=P3(i_im_jp); vip_jp=P3(i_ip_jp)
      CASE (4); vim_jm=P4(i_im_jm); vip_jm=P4(i_ip_jm); vim_jp=P4(i_im_jp); vip_jp=P4(i_ip_jp)
      CASE (5); vim_jm=P5(i_im_jm); vip_jm=P5(i_ip_jm); vim_jp=P5(i_im_jp); vip_jp=P5(i_ip_jp)
      CASE (6); vim_jm=P6(i_im_jm); vip_jm=P6(i_ip_jm); vim_jp=P6(i_im_jp); vip_jp=P6(i_ip_jp)
      CASE DEFAULT
                vim_jm=P7(i_im_jm); vip_jm=P7(i_ip_jm); vim_jp=P7(i_im_jp); vip_jp=P7(i_ip_jp)
    END SELECT
  END SUBROUTINE stark_profile_interp


! ============================================================================
!  hydrogen_line_profile -- hydrogen line profile including Stark broadening and fine structure.
!
!  Computes the opacity profile for a single hydrogen line transition N->M
!  at depth point J, for a wavelength offset DELW from line centre.
!
!  Arguments:
!    N      -- lower principal quantum number
!    M      -- upper principal quantum number
!    J      -- atmospheric depth index
!    DELW   -- wavelength offset from line centre (nm, REAL*8)
!    DOPPH  -- Doppler half-width array (per depth; REAL*4, dimension kw)
!
!  Returns: opacity profile value (REAL*4), normalised as for voigt_profile.
!
!  Physics (from Deane Peterson / Kurucz):
!    Doppler core: fine-structure components convolved with Gaussian
!    Lorentz wings: resonance, radiation, van der Waals broadening
!    Lyman alpha: special treatment with H2/H2+ quasi-molecular cutoffs
!                 (Allard & Koester 1992; Allard et al. 1998)
!    Stark wings: quasistatic ion field profile S(beta,p) via stark_quasistatic_profile
! ============================================================================
  FUNCTION hydrogen_line_profile(n, m, j, delw, dopph) RESULT(result)
    INTEGER,  INTENT(IN) :: n, m, j
    REAL(8),  INTENT(IN) :: delw
    REAL(4),  INTENT(IN) :: dopph(kw)
    REAL(4)              :: result

    ! Lyman alpha H2+ quasi-molecular cutoff table
    ! Delta wavenumber = -15000 + 100*(i-1), i=1..111, up to -4000 cm^{-1}
    REAL(4), PARAMETER :: CUTOFFH2PLUS(111) = [ &
      -15.14,-15.06,-14.97,-14.88,-14.80,-14.71,-14.62,-14.53, &
      -14.44,-14.36,-14.27,-14.18,-14.09,-14.01,-13.92,-13.83, &
      -13.74,-13.65,-13.57,-13.48,-13.39,-13.30,-13.21,-13.13, &
      -13.04,-12.95,-12.86,-12.77,-12.69,-12.60,-12.51,-12.40, &
      -12.29,-12.15,-12.02,-11.90,-11.76,-11.63,-11.53,-11.41, &
      -11.30,-11.22,-11.15,-11.09,-11.07,-11.06,-11.07,-11.09, &
      -11.12,-11.16,-11.19,-11.21,-11.24,-11.27,-11.30,-11.33, &
      -11.36,-11.39,-11.42,-11.45,-11.48,-11.48,-11.48,-11.48, &
      -11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48, &
      -11.48,-11.48,-11.48,-11.48,-11.41,-11.40,-11.39,-11.38, &
      -11.37,-11.36,-11.35,-11.34,-11.33,-11.32,-11.30,-11.29, &
      -11.28,-11.27,-11.27,-11.27,-11.26,-11.25,-11.24,-11.23, &
      -11.22,-11.21,-11.20,-11.19,-11.18,-11.17,-11.15,-11.14, &
      -11.13,-11.12,-11.11,-11.10,-11.09,-11.08,-11.07 ]

    ! Lyman alpha H2 quasi-molecular cutoff table
    ! Delta wavenumber = -22000 + 200*(i-1), i=1..91, up to -4000 cm^{-1}
    REAL(4), PARAMETER :: CUTOFFH2(91) = [ &
      -13.64,-13.52,-13.39,-13.27,-13.14,-13.01,-12.87,-12.74, &
      -12.63,-12.56,-12.51,-12.48,-12.47,-12.49,-12.52,-12.55, &
      -12.57,-12.61,-12.65,-12.69,-12.72,-12.76,-12.79,-12.82, &
      -12.84,-12.85,-12.87,-12.90,-12.93,-12.94,-12.93,-12.95, &
      -12.95,-12.96,-12.97,-12.96,-12.96,-12.95,-12.95,-12.96, &
      -12.98,-12.99,-12.95,-12.96,-13.00,-13.00,-12.98,-12.97, &
      -13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, &
      -13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, &
      -13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-12.89,-12.88, &
      -12.87,-12.86,-12.85,-12.84,-12.83,-12.81,-12.80,-12.79, &
      -12.78,-12.76,-12.74,-12.72,-12.70,-12.68,-12.65,-12.62, &
      -12.59,-12.56,-12.53 ]

    ! Einstein A coefficients (sum over upper and lower levels)
    REAL(4), PARAMETER :: ASUMLYMAN(100) = [ &
       0.000E+00, 6.265E+08, 1.897E+08, 8.126E+07, 4.203E+07, 2.450E+07, &
       1.236E+07, 8.249E+06, 5.782E+06, 4.208E+06, 3.158E+06, 2.430E+06, &
       1.910E+06, 1.567E+06, 1.274E+06, 1.050E+06, 8.752E+05, 7.373E+05, &
       6.269E+05, 5.375E+05, 4.643E+05, 4.038E+05, 3.534E+05, 3.111E+05, &
       2.752E+05, 2.447E+05, 2.185E+05, 1.959E+05, 1.763E+05, 1.593E+05, &
       1.443E+05, 1.312E+05, 1.197E+05, 1.094E+05, 1.003E+05, 9.216E+04, &
       8.489E+04, 7.836E+04, 7.249E+04, 6.719E+04, 6.239E+04, 5.804E+04, &
       5.408E+04, 5.048E+04, 4.719E+04, 4.418E+04, 4.142E+04, 3.888E+04, &
       3.655E+04, 3.440E+04, 3.242E+04, 3.058E+04, 2.888E+04, 2.731E+04, &
       2.585E+04, 2.449E+04, 2.322E+04, 2.204E+04, 2.094E+04, 1.991E+04, &
       1.894E+04, 1.804E+04, 1.720E+04, 1.640E+04, 1.566E+04, 1.496E+04, &
       1.430E+04, 1.368E+04, 1.309E+04, 1.254E+04, 1.201E+04, 1.152E+04, &
       1.105E+04, 1.061E+04, 1.019E+04, 9.796E+03, 9.419E+03, 9.061E+03, &
       8.721E+03, 8.398E+03, 8.091E+03, 7.799E+03, 7.520E+03, 7.255E+03, &
       7.002E+03, 6.760E+03, 6.530E+03, 6.310E+03, 6.100E+03, 5.898E+03, &
       5.706E+03, 5.522E+03, 5.346E+03, 5.177E+03, 5.015E+03, 4.860E+03, &
       4.711E+03, 4.569E+03, 4.432E+03, 4.300E+03 ]

    REAL(4), PARAMETER :: ASUM(100) = [ &
       0.000E+00, 4.696E+08, 9.980E+07, 3.017E+07, 1.155E+07, 5.189E+06, &
       2.616E+06, 1.437E+06, 8.444E+05, 5.234E+05, 3.389E+05, 2.275E+05, &
       1.575E+05, 1.120E+05, 8.142E+04, 6.040E+04, 4.560E+04, 3.496E+04, &
       2.719E+04, 2.141E+04, 1.711E+04, 1.377E+04, 1.119E+04, 9.166E+03, &
       7.572E+03, 6.341E+03, 5.338E+03, 4.523E+03, 3.854E+03, 3.302E+03, &
       2.844E+03, 2.460E+03, 2.138E+03, 1.866E+03, 1.635E+03, 1.438E+03, &
       1.269E+03, 1.124E+03, 9.983E+02, 8.894E+02, 7.947E+02, 7.120E+02, &
       6.396E+02, 5.759E+02, 5.198E+02, 4.703E+02, 4.263E+02, 3.873E+02, &
       3.526E+02, 3.215E+02, 2.938E+02, 2.689E+02, 2.465E+02, 2.264E+02, &
       2.082E+02, 1.918E+02, 1.769E+02, 1.634E+02, 1.512E+02, 1.400E+02, &
       1.298E+02, 1.206E+02, 1.121E+02, 1.043E+02, 9.720E+01, 9.066E+01, &
       8.465E+01, 7.912E+01, 7.403E+01, 6.933E+01, 6.498E+01, 6.097E+01, &
       5.725E+01, 5.381E+01, 5.061E+01, 4.765E+01, 4.489E+01, 4.232E+01, &
       3.994E+01, 3.771E+01, 3.563E+01, 3.369E+01, 3.188E+01, 3.019E+01, &
       2.860E+01, 2.712E+01, 2.572E+01, 2.442E+01, 2.319E+01, 2.204E+01, &
       2.096E+01, 1.994E+01, 1.898E+01, 1.808E+01, 1.722E+01, 1.642E+01, &
       1.566E+01, 1.495E+01, 1.427E+01, 1.363E+01 ]

    ! Fine-structure component offsets and weights (Kurucz)
    REAL(4), PARAMETER :: STALPH(34) = [ &
      -730., 370., 188., 515., 327., 619.,-772.,-473.,-369., 120., &
       256., 162., 285.,-161., -38.3,  6.82,-174.,-147.,-101., -77.5, &
        55., 126.,  75., 139., -60.,   3.7,  27., -69., -42., -18., &
        -5.5,  -9.1, -33., -24. ]
    REAL(4), PARAMETER :: STWTAL(34) = [ &
      1.,2.,1.,2.,1.,2.,1.,2.,3.,1.,2.,1.,2.,1.,4.,6.,1.,2., &
      3.,4.,1.,2.,1.,2.,1.,4.,6.,1.,7.,6.,4.,4.,4.,5. ]
    INTEGER, PARAMETER :: ISTAL(4) = [ 1, 3, 10, 21 ]
    INTEGER, PARAMETER :: LNGHAL(4) = [ 2, 7, 11, 14 ]

    REAL(4), PARAMETER :: STCOMP(5,4) = RESHAPE( [ &
       0.,  0.,   0.,   0., 0., &
     468.,576.,-522.,   0., 0., &
     260.,290., -33.,-140., 0., &
     140.,150.,  18., -27.,-51. ], [5,4] )
    REAL(4), PARAMETER :: STCPWT(5,4) = RESHAPE( [ &
      1.,0.,0.,0.,0., &
      1.,1.,2.,0.,0., &
      1.,1.,4.,3.,0., &
      1.,1.,4.,6.,4. ], [5,4] )
    INTEGER, PARAMETER :: LNCOMP(4) = [ 1, 3, 4, 5 ]

    ! XKNM table (Kurucz) for principal quantum numbers n<=4, delta_m<=3
    REAL(4), PARAMETER :: XKNMTB(4,3) = RESHAPE( [ &
      0.0001716, 0.009019, 0.1001, 0.5820, &
      0.0005235, 0.01772,  0.171,  0.866,  &
      0.0008912, 0.02507,  0.223,  1.02    ], [4,3] )

    INTEGER, SAVE :: itemp1 = 0, n1 = 0, m1 = 0
    INTEGER       :: k, i, ifins, ipos, nwid, ifcore, icut
    INTEGER       :: mmn
    REAL(4)       :: xn, xm, xn2, xm2, xmn2, xm2mn2, gnm, xknm
    REAL(4)       :: y1num, y1wht, freqnm, dbeta, wavenm
    REAL(4)       :: c1con, c2con, radamp, resont, vdw, hwvdw
    REAL(4)       :: stark, wl4, freq, del
    REAL(4)       :: hwstk, hwlor, hwdop, hwres, hwrad, hfwid, dop
    REAL(4)       :: finest(14), finswt(14)
    REAL(4)       :: hprofres, hprofvdw, hprofrad, hproflor, hhw, top
    REAL(4)       :: wty1, y1scal, c1_loc, c2_loc, g1, gnot
    REAL(4)       :: beta_val, y1, y2, gam, prqs, f, p1, fns
    REAL(4)       :: cutoff_val, cutfreq, spacing
    REAL(4)       :: freq22000, freq15000, beta4000, prqsp4000, cutoff4000
    REAL(4)       :: d, hh1, hh2, hh3, hh4
    ! Per-depth pre-computed arrays (saved, recomputed when ITEMP changes)
    REAL(4), SAVE :: pp_d(kw), fo(kw), gcon1(kw), gcon2(kw)
    REAL(4), SAVE :: y1b(kw), y1s(kw), c1d(kw), c2d(kw)
    REAL(4), SAVE :: t3nhe(kw), t3nh2(kw)
    REAL(4), SAVE :: radamp_s, resont_s, vdw_s, stark_s
    REAL(4), SAVE :: c1con_s, c2con_s, dbeta_s, wavenm_s, freqnm_s
    REAL(4), SAVE :: y1num_s, y1wht_s
    REAL(4), SAVE :: hwvdw_s, hwrad_s
    INTEGER, SAVE :: ifins_s, nwid_s
    REAL(4), SAVE :: finest_s(14), finswt_s(14)
    REAL(4)       :: xne16, t4, t43

    REAL(4), PARAMETER :: RYDH = 3.2880515E15

    ! --- Re-compute depth-dependent quantities if temperature index changed ---
    IF (itemp /= itemp1) THEN
      itemp1 = itemp
      DO k = 1, NRHOX
        xne16     = XNE(k)**0.1666667
        pp_d(k)   = xne16 * 0.08989 / SQRT(T(k))
        fo(k)     = xne16**4 * 1.25E-9
        y1b(k)    = 2.0 / (1.0 + 0.012/T(k) * SQRT(XNE(k)/T(k)))
        t4        = T(k) / 10000.0
        t43       = t4**0.3
        y1s(k)    = t43 / xne16
        t3nhe(k)  = t43 * XNFHE(k,1)
        t3nh2(k)  = t43 * XNFH2(k)
        c1d(k)    = fo(k) * 78940.0 / T(k)
        c2d(k)    = fo(k)**2 / 5.96E-23 / XNE(k)
        gcon1(k)  = 0.2 + 0.09*SQRT(t4) / (1.0 + XNE(k)/1.E13)
        gcon2(k)  = 0.2 / (1.0 + XNE(k)/1.E15)
      END DO
    END IF

    ! --- Re-compute line-dependent quantities if N or M changed ---
    IF (n /= n1 .OR. m /= m1) THEN
      n1  = n;  m1 = m
      mmn = m - n
      xn  = REAL(n);  xm = REAL(m)
      xn2 = xn*xn;    xm2 = xm*xm
      xmn2   = xm2 * xn2
      xm2mn2 = xm2 - xn2
      gnm  = xm2mn2 / xmn2

      IF (mmn <= 3 .AND. n <= 4) THEN
        xknm = XKNMTB(n, mmn)
      ELSE
        xknm = 5.5E-5 / gnm * xmn2 / (1.0 + 0.13/REAL(mmn))
      END IF

      y1num_s = 320.0
      IF (m == 2) y1num_s = 550.0
      IF (m == 3) y1num_s = 380.0
      y1wht_s = 1.E13
      IF (mmn <= 3) y1wht_s = 1.E14
      IF (mmn <= 2 .AND. n <= 2) THEN
        ! Y1WTM table: (1->2)=1e18, (1->3)=1e17, (2->2)=1e16, (2->3)=1e14
        IF      (n==1 .AND. mmn==1) THEN; y1wht_s = 1.E18
        ELSE IF (n==1 .AND. mmn==2) THEN; y1wht_s = 1.E17
        ELSE IF (n==2 .AND. mmn==1) THEN; y1wht_s = 1.E16
        ELSE IF (n==2 .AND. mmn==2) THEN; y1wht_s = 1.E14
        END IF
      END IF

      freqnm_s = RYDH * gnm
      dbeta_s  = REAL(CLIGHT_AHZ,4) / freqnm_s**2 / xknm
      wavenm_s = REAL(CLIGHT_AHZ,4) / freqnm_s
      c1con_s  = xknm / wavenm_s * gnm * xm2mn2
      c2con_s  = (xknm / wavenm_s)**2

      radamp_s = ASUM(n) + ASUM(m)
      IF (n == 1) radamp_s = ASUMLYMAN(m)
      radamp_s = radamp_s / 12.5664 / freqnm_s

      resont_s = hydrogen_oscillator_strength(1,m)/xm/(1.0-1.0/xm2)
      IF (n /= 1) resont_s = resont_s + hydrogen_oscillator_strength(1,n)/xn/(1.0-1.0/xn2)
      resont_s = resont_s * 3.579E-24 / gnm

      vdw_s    = 4.45E-26 / gnm * (xm2*(7.0*xm2 + 5.0))**0.4
      stark_s  = 1.6678E-18 * freqnm_s * xknm

      ! Fine structure
      IF (n > 4 .OR. m > 10) THEN
        ifins_s    = 1
        finest_s(1) = 0.0
        finswt_s(1) = 1.0
      ELSE IF (mmn == 1) THEN
        ! Alpha line: explicit fine structure
        ifins_s = LNGHAL(n)
        ipos    = ISTAL(n)
        DO i = 1, ifins_s
          finest_s(i) = STALPH(ipos-1+i) * 1.E7
          finswt_s(i) = STWTAL(ipos-1+i) / xn2 / 3.0
        END DO
      ELSE
        ! Non-alpha: M=inf pattern
        ifins_s = LNCOMP(n)
        DO i = 1, ifins_s
          finest_s(i) = STCOMP(i,n) * 1.E7
          finswt_s(i) = STCPWT(i,n) / xn2
        END DO
      END IF
    END IF

    ! --- Now compute profile at this (j, delw) ---
    ! Recover saved line quantities into locals for readability
    freqnm = freqnm_s;  dbeta  = dbeta_s;  wavenm = wavenm_s
    c1con  = c1con_s;   c2con  = c2con_s;  radamp = radamp_s
    resont = resont_s;  vdw    = vdw_s;    stark  = stark_s
    y1num  = y1num_s;   y1wht  = y1wht_s
    ifins  = ifins_s
    finest(1:ifins) = finest_s(1:ifins)
    finswt(1:ifins) = finswt_s(1:ifins)

    ! Wavelength and frequency at this offset (REAL*8 delw, converted)
    wl4  = REAL(wavenm + delw*10.0)    ! wavelength in Angstrom, then /10 -> nm below
    freq = REAL(CLIGHT_AHZ,4) / wl4
    del  = ABS(freq - freqnm)

    ! Broadening half-widths (dimensionless nu/nu_0)
    hwstk = stark * fo(j)
    hwvdw = vdw * t3nhe(j) + 2.0*vdw * t3nh2(j)
    hwrad = radamp
    hwres = resont * XNFPH(j,1) * 2.0
    hwlor = hwres + hwvdw + hwrad
    hwdop = dopph(j)

    ! Identify dominant broadening mechanism
    nwid = 1
    IF (hwdop < hwstk .OR. hwdop < hwlor) THEN
      nwid = 2
      IF (hwlor < hwstk) nwid = 3
    END IF
    hfwid = freqnm * MAX(hwdop, hwlor, hwstk)

    result  = 0.0
    ifcore  = 0
    IF (del <= hfwid) ifcore = 1
    dop = freqnm * hwdop

    ! ---- DOPPLER CORE ----
    IF (ifcore == 1 .AND. nwid == 1) THEN
      DO i = 1, ifins
        d = ABS(freq - freqnm - finest(i)) / dop
        IF (d <= 7.0) result = result + EXP(-d*d) * finswt(i)
      END DO
      RETURN
    END IF

    ! ---- LORENTZ WINGS ----
    IF (ifcore == 1 .AND. nwid == 2 .OR. nwid == 2 .AND. ifcore == 0 .OR. &
        nwid == 3 .AND. ifcore == 0) THEN
      ! Lyman alpha special treatment
      IF (n == 1 .AND. m == 2) THEN
        ASSOCIATE(hwres4 => hwres * 4.0)
          hwlor    = hwres4 + hwvdw + hwrad
          hhw      = freqnm * hwlor
          IF (freq > (82259.10 - 4000.0) * REAL(CLIGHT_CMS,4)) THEN
            hprofres = hwres4 * freqnm / 3.14159 / (del**2 + hhw**2) * 1.77245 * dop
          ELSE
            ! Far red wing: H2 quasi-molecular cutoff (Allard & Koester 1992)
            cutoff_val = 0.0
            spacing    = 200.0 * REAL(CLIGHT_CMS,4)
            freq22000  = (82259.10 - 22000.0) * REAL(CLIGHT_CMS,4)
            IF (freq >= 50000.0 * REAL(CLIGHT_CMS,4)) THEN
              IF (freq < freq22000) THEN
                cutoff_val = (CUTOFFH2(2)-CUTOFFH2(1))/spacing*(freq-freq22000) + CUTOFFH2(1)
              ELSE
                icut = INT((freq - freq22000) / spacing)
                cutfreq = icut*spacing + freq22000
                cutoff_val = (CUTOFFH2(icut+2)-CUTOFFH2(icut+1))/spacing*(freq-cutfreq) + CUTOFFH2(icut+1)
              END IF
              cutoff_val = (10.0**(cutoff_val-14.0)) * XNFPH(j,1) * 2.0 / REAL(CLIGHT_CMS,4)
            END IF
            hprofres = cutoff_val * 1.77245 * dop
          END IF
        END ASSOCIATE
        ! Radiation damping contribution
        hprofrad = 0.0
        IF (freq > 2.4190611E15 .AND. freq < 0.77*3.288051E15) &
          hprofrad = hwrad * freqnm / 3.14159 / (del**2 + hhw**2) * 1.77245 * dop
        ! Van der Waals
        hprofvdw = hwvdw * freqnm / 3.14159 / (del**2 + hhw**2) * 1.77245 * dop
        IF (freq < 1.8E15) hprofvdw = 0.0
        hproflor = hprofres + hprofrad + hprofvdw
        result   = result + hproflor
        IF (ifcore == 1) RETURN
      ELSE
        ! Standard Lorentz wing
        hhw = freqnm * hwlor
        top = hhw
        ! Lyman beta/gamma/delta: suppress radiation near series limit
        IF (n == 1 .AND. m <= 5) THEN
          IF (m==3 .AND. freq > 0.885*3.288051E15 .AND. freq < 0.890*3.288051E15) &
            top = hhw - freqnm*hwrad
          IF (m==4 .AND. freq > 0.936*3.288051E15 .AND. freq < 0.938*3.288051E15) &
            top = hhw - freqnm*hwrad
          IF (m==5 .AND. freq > 0.959*3.288051E15 .AND. freq < 0.961*3.288051E15) &
            top = hhw - freqnm*hwrad
        END IF
        hproflor = top / 3.14159 / (del**2 + hhw**2) * 1.77245 * dop
        result   = result + hproflor
        IF (ifcore == 1) RETURN
      END IF
    END IF

    ! ---- STARK WINGS ----
    wty1    = 1.0 / (1.0 + XNE(j) / y1wht)
    y1scal  = y1num * y1s(j) * wty1 + y1b(j) * (1.0 - wty1)
    c1_loc  = c1d(j) * c1con * y1scal
    c2_loc  = c2d(j) * c2con
    g1      = 6.77 * SQRT(c1_loc)
    gnot    = g1 * MAX(0.0, 0.2114 + LOG(SQRT(c2_loc)/c1_loc)) &
                 * (1.0 - gcon1(j) - gcon2(j))
    beta_val = ABS(del) / fo(j) * dbeta
    y1      = c1_loc * beta_val
    y2      = c2_loc * beta_val**2
    gam     = gnot
    IF (.NOT. (y2 <= 1.E-4 .AND. y1 <= 1.E-5)) THEN
      gam = g1 * (0.5*EXP(-MIN(80.0, y1)) + expint1(y1) - 0.5*expint1(y2)) &
              * (1.0 - gcon1(j)/(1.0+(90.0*y1)**3) - gcon2(j)/(1.0+2000.0*y1))
      IF (gam <= 1.E-20) gam = 0.0
    END IF

    prqs = stark_quasistatic_profile(beta_val, pp_d(j), n, m)

    ! Lyman alpha: split quasistatic profile, add H2+ cutoff
    IF (m == 2) THEN
      prqs = prqs * 0.5
      cutoff_val = 0.0
      IF (freq >= (82259.10-20000.0)*REAL(CLIGHT_CMS,4) .AND. &
          freq <= (82259.10- 4000.0)*REAL(CLIGHT_CMS,4)) THEN
        spacing   = 100.0 * REAL(CLIGHT_CMS,4)
        freq15000 = (82259.10 - 15000.0) * REAL(CLIGHT_CMS,4)
        IF (freq < freq15000) THEN
          cutoff_val = (CUTOFFH2PLUS(2)-CUTOFFH2PLUS(1))/spacing*(freq-freq15000) + CUTOFFH2PLUS(1)
        ELSE
          icut     = INT((freq - freq15000) / spacing)
          cutfreq  = icut*spacing + freq15000
          cutoff_val = (CUTOFFH2PLUS(icut+2)-CUTOFFH2PLUS(icut+1))/spacing*(freq-cutfreq) + CUTOFFH2PLUS(icut+1)
        END IF
        cutoff_val = (10.0**(cutoff_val-14.0)) / REAL(CLIGHT_CMS,4) * XNFPH(j,2)
        result = result + cutoff_val * 1.77245 * dop
      ELSE IF (freq > (82259.10-4000.0)*REAL(CLIGHT_CMS,4)) THEN
        beta4000    = 4000.0*REAL(CLIGHT_CMS,4) / fo(j) * dbeta
        prqsp4000   = stark_quasistatic_profile(beta4000, pp_d(j), n, m) * 0.5 / fo(j) * dbeta
        cutoff4000  = (10.0**(-11.07-14.0)) / REAL(CLIGHT_CMS,4) * XNFPH(j,2)
        result = result + cutoff4000/prqsp4000 * prqs/fo(j) * dbeta * 1.77245 * dop
      END IF
    END IF

    f = 0.0
    IF (gam > 0.0) f = gam / 3.14159 / (gam**2 + beta_val**2)
    p1  = (0.9*y1)**2
    fns = (p1 + 0.03*SQRT(y1)) / (p1 + 1.0)
    result = result + (prqs*(1.0 + fns) + f) / fo(j) * dbeta * 1.77245 * dop

  END FUNCTION hydrogen_line_profile



! ============================================================================
!  he1_generic_profile -- dispatcher for He I line profiles.
!
!  Selects the appropriate detailed broadening function based on line
!  wavelength (encoded as an integer nm value).  Falls back to he1_griem_profile
!  for lines not explicitly handled.
!
!  Arguments:
!    J       -- atmospheric depth index
!    WAVE    -- current wavelength (nm)
!    WL      -- line-centre wavelength (nm)
!    DOPWL   -- Doppler width in wavelength units (nm)
!    GAMMAR  -- radiation damping constant
!    GAMMAS  -- Stark damping constant
! ============================================================================
  FUNCTION he1_generic_profile(j, wave, wl, dopwl, gammar, gammas) RESULT(result)
    INTEGER, INTENT(IN) :: j
    REAL(4), INTENT(IN) :: wave, wl, dopwl, gammar, gammas
    REAL(4)             :: result
    INTEGER             :: line

    line = INT(wl + 1.0)

    SELECT CASE (line)
      CASE (448)
        result = he1_4471_profile(j, wave, wl, dopwl)
        RETURN
      CASE (403)
        ! Distinguish 402.6 from 402.3 nm
        IF (INT(wl + 0.4) == 402) THEN
          ! 402.3 nm: falls through to he1_griem_profile check below
        ELSE
          result = he1_4026_profile(j, wave, wl, dopwl)
          RETURN
        END IF
      CASE (439)
        result = he1_4387_profile(j, wave, wl, dopwl)
        RETURN
      CASE (493)
        result = he1_4921_profile(j, wave, wl, dopwl)
        RETURN
    END SELECT

    ! Other listed He I lines with he1_dimitri_profile broadening
    IF (line==382 .OR. line==387 .OR. line==393 .OR. &
        line==401 .OR. line==403 .OR. line==415 .OR. line==417) THEN
      result = he1_dimitri_profile(j, wave, wl, dopwl, gammar, gammas)
      RETURN
    END IF

    ! Default: he1_griem_profile broadening
    result = he1_griem_profile(j, wave, wl, dopwl, gammar, gammas)

  END FUNCTION he1_generic_profile


! ============================================================================
!  he1_4471_profile -- He I 447.1 nm line profile.
!
!  Uses isolated-line Stark broadening parameters from:
!    Barnard, Cooper & Smith, J.Q.S.R.T. 14, 1025 (1974).
!  Fine-structure splitting included explicitly.
!  Falls back to read_he1_stark_tables tabulated profile for high electron densities.
! ============================================================================
  FUNCTION he1_4471_profile(j, wave, wl, dopwl) RESULT(result)
    INTEGER, INTENT(IN) :: j
    REAL(4), INTENT(IN) :: wave, wl, dopwl
    REAL(4)             :: result

    ! Broadening parameters at Ne = 1e13 cm^{-3} (Barnard et al. 1974)
    REAL(4), PARAMETER :: WS(4)   = [ 0.001460, 0.001269, 0.001079, 0.000898 ]
    REAL(4), PARAMETER :: DS(4)   = [ 0.036, -0.005, -0.026, -0.034 ]
    REAL(4), PARAMETER :: ALFS(4) = [ 0.107,  0.119,  0.134,  0.154 ]
    REAL(4), PARAMETER :: FORB1(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: FORB2(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: TS(4)   = [ 5.0E3, 1.0E4, 2.0E4, 4.0E4 ]
    REAL(4), PARAMETER :: DEN     = 1.0E13
    ! Fine structure splittings (nm)
    REAL(4), PARAMETER :: DLF1   = -0.150
    REAL(4), PARAMETER :: DLF2   =  0.0
    ! Fine structure positions for the 2p3P -> 4d3D blend (nm offsets from WL)
    REAL(4), PARAMETER :: FS1 = -0.0184, FS2 = 0.0013, FS3 = 0.0010, &
                           FS4 = -0.0029, FS5 = -0.0025
    REAL(4), PARAMETER :: FW1 = 1.0/9.0, FW2 = 1.0/12.0, FW3 = 0.25, &
                           FW4 = 1.0/180.0, FW5 = 11.0/20.0
    REAL(4) :: phihe, wtot, dtot, wwd, a_damp, v, sentinel

    result   = 0.0
    sentinel = 0.0
    CALL he1_isolated_profile(j, wave, wl, dopwl, WS, DS, ALFS, FORB1, FORB2, &
                             TS, DEN, DLF1, DLF2, 1.0E13, 1, sentinel)
    IF (sentinel < 0.0) THEN
      CALL read_he1_stark_tables(1, j, T(j), XNFPH(j,2), XNFHE(j,2), XNE(j), wave-wl, phihe)
      result = 1.772453 * phihe * dopwl * 10.0
      RETURN
    END IF
    ! Apply fine structure splitting
    CALL he1_stark_widths(j, wave, wl, WS, DS, ALFS, TS, DEN, wtot, dtot, a_damp)
    wwd    = wave - wl - dtot
    result = voigt_profile(ABS(wwd - FS1)/dopwl, a_damp) * FW1 + &
             voigt_profile(ABS(wwd + FS2)/dopwl, a_damp) * FW2 + &
             voigt_profile(ABS(wwd + FS3)/dopwl, a_damp) * FW3 + &
             voigt_profile(ABS(wwd + FS4)/dopwl, a_damp) * FW4 + &
             voigt_profile(ABS(wwd + FS5)/dopwl, a_damp) * FW5
    IF (FORB1(1) > 0.0) THEN
      v      = ABS((wave - wl) - DLF1) / dopwl
      result = result + FORB1(1) * voigt_profile(v, a_damp)  ! FORB1 is 0 per DATA
    END IF

  END FUNCTION he1_4471_profile


! ============================================================================
!  Shared helper: compute Stark width and shift for isolated He I lines.
!
!  Returns result < 0 as a sentinel when XNE(J) exceeds the density
!  threshold DENHI (meaning read_he1_stark_tables should be used instead).
! ============================================================================
  SUBROUTINE he1_isolated_profile(j, wave, wl, dopwl, ws, ds, alfs, forb1, forb2, &
                                  ts, den, dlf1, dlf2, denhi, line_flag, result)
    INTEGER, INTENT(IN)  :: j, line_flag
    REAL(4), INTENT(IN)  :: wave, wl, dopwl
    REAL(4), INTENT(IN)  :: ws(4), ds(4), alfs(4), forb1(4), forb2(4)
    REAL(4), INTENT(IN)  :: ts(4), den, dlf1, dlf2, denhi
    REAL(4), INTENT(OUT) :: result

    REAL(4) :: e, temp, x, xx, w_damp, d_ratio, alf
    REAL(4) :: f1, f2, xnfhp, xnfhep, vm1, rhom, sigma, wtot, dtot, a_damp, v
    INTEGER :: it

    e      = XNE(j)
    temp   = T(j)
    xnfhp  = XNFPH(j,2)
    xnfhep = XNFHE(j,2)
    temp   = MAX(temp, 5.0E3)
    temp   = MIN(temp, 4.0E4)

    IF (e > denhi) THEN
      result = -1.0   ! sentinel: use read_he1_stark_tables
      RETURN
    END IF

    ! Temperature interpolation index
    it = 2
    DO WHILE (it < 4 .AND. ts(it) <= temp)
      it = it + 1
    END DO
    x  = (temp - ts(it-1)) / (ts(it) - ts(it-1))
    xx = e / den

    w_damp  = xx * (x*ws(it) + (1.0-x)*ws(it-1))
    d_ratio = x*ds(it) + (1.0-x)*ds(it-1)
    alf     = xx**0.25 * (x*alfs(it) + (1.0-x)*alfs(it-1))
    f1      = xx * (x*forb1(it) + (1.0-x)*forb1(it-1))
    f2      = xx * (x*forb2(it) + (1.0-x)*forb2(it-1))

    xx     = xnfhp / e
    vm1    = 8.78 * (xx + 2.0*(1.0-xx)) / SQRT(temp)
    rhom   = 1.0 / (4.19*e)**(1.0/3.0)
    sigma  = 1.885E14 * w_damp * rhom * vm1 / (wl*10.0)**2
    x      = alf**(8.0/9.0) / sigma**(1.0/3.0)
    wtot   = w_damp * (1.0 + 1.36*x) * 0.1    ! Angstrom -> nm
    dtot   = w_damp * d_ratio * (1.0 + 2.36*x/ABS(d_ratio)) * 0.1

    a_damp = wtot / dopwl
    ! Caller will apply fine structure; just return a sentinel >= 0
    result = a_damp  ! positive: means "use fine-structure path"

  END SUBROUTINE he1_isolated_profile

  ! Recompute widths for fine-structure application
  SUBROUTINE he1_stark_widths(j, wave, wl, ws, ds, alfs, ts, den, wtot, dtot, a_damp)
    INTEGER, INTENT(IN)  :: j
    REAL(4), INTENT(IN)  :: wave, wl, ws(4), ds(4), alfs(4), ts(4), den
    REAL(4), INTENT(OUT) :: wtot, dtot, a_damp
    REAL(4) :: e, temp, x, xx, w_d, d_r, alf, xnfhp, vm1, rhom, sigma, dopwl_loc
    INTEGER :: it
    e      = XNE(j);  temp = MAX(MIN(T(j), 4.0E4), 5.0E3)
    xnfhp  = XNFPH(j,2)
    it = 2
    DO WHILE (it < 4 .AND. ts(it) <= temp)
      it = it + 1
    END DO
    x     = (temp - ts(it-1)) / (ts(it) - ts(it-1));  xx = e/den
    w_d   = xx*(x*ws(it)+(1.0-x)*ws(it-1))
    d_r   = x*ds(it)+(1.0-x)*ds(it-1)
    alf   = xx**0.25*(x*alfs(it)+(1.0-x)*alfs(it-1))
    xx    = xnfhp/e
    vm1   = 8.78*(xx+2.0*(1.0-xx))/SQRT(temp)
    rhom  = 1.0/(4.19*e)**(1.0/3.0)
    sigma = 1.885E14*w_d*rhom*vm1/(wl*10.0)**2
    x     = alf**(8.0/9.0)/sigma**(1.0/3.0)
    wtot  = w_d*(1.0+1.36*x)*0.1
    dtot  = w_d*d_r*(1.0+2.36*x/ABS(d_r))*0.1
    dopwl_loc = REAL(DOPPLE(7))*wl
    a_damp = wtot/dopwl_loc
  END SUBROUTINE he1_stark_widths


! ============================================================================
!  he1_4026_profile, he1_4387_profile, he1_4921_profile -- He I line profiles at 402.6, 438.7, 492.1 nm.
!
!  Same structure as he1_4471_profile; parameters from he1_griem_profile (1974) and
!  Barnard, Cooper & Smith (1975).  Fine structure differs per line.
! ============================================================================
  FUNCTION he1_4026_profile(j, wave, wl, dopwl) RESULT(result)
    INTEGER, INTENT(IN) :: j
    REAL(4), INTENT(IN) :: wave, wl, dopwl
    REAL(4)             :: result
    REAL(4), PARAMETER :: WS(4)   = [ 4.04,  3.49,  2.96,  2.47 ]
    REAL(4), PARAMETER :: DS(4)   = [ 0.1339, 0.0960, 0.0780, 0.0709 ]
    REAL(4), PARAMETER :: ALFS(4) = [ 0.969,  1.083,  1.225,  1.403 ]
    REAL(4), PARAMETER :: FORB1(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: FORB2(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: TS(4)   = [ 5.0E3, 1.0E4, 2.0E4, 4.0E4 ]
    REAL(4), PARAMETER :: DEN = 1.0E16, DLF1 = -0.080, DLF2 = 0.0
    REAL(4), PARAMETER :: FS1=-0.0148, FS2=0.0012, FS3=0.0011, FS4=-0.0025, FS5=-0.0023
    REAL(4), PARAMETER :: FW1=1.0/9.0, FW2=1.0/12.0, FW3=0.25, FW4=1.0/180.0, FW5=11.0/20.0
    REAL(4) :: phihe, wtot, dtot, a_damp, wwd, sentinel

    result = 0.0
    CALL he1_isolated_profile(j, wave, wl, dopwl, WS, DS, ALFS, FORB1, FORB2, &
                             TS, DEN, DLF1, DLF2, 1.0E14, 2, sentinel)
    IF (sentinel < 0.0) THEN
      CALL read_he1_stark_tables(2, j, T(j), XNFPH(j,2), XNFHE(j,2), XNE(j), wave-wl, phihe)
      result = 1.772453 * phihe * dopwl * 10.0
      RETURN
    END IF
    CALL he1_stark_widths(j, wave, wl, WS, DS, ALFS, TS, DEN, wtot, dtot, a_damp)
    wwd    = wave - wl - dtot
    result = voigt_profile(ABS(wwd-FS1)/dopwl, a_damp)*FW1 + &
             voigt_profile(ABS(wwd+FS2)/dopwl, a_damp)*FW2 + &
             voigt_profile(ABS(wwd+FS3)/dopwl, a_damp)*FW3 + &
             voigt_profile(ABS(wwd+FS4)/dopwl, a_damp)*FW4 + &
             voigt_profile(ABS(wwd+FS5)/dopwl, a_damp)*FW5
  END FUNCTION he1_4026_profile

  FUNCTION he1_4387_profile(j, wave, wl, dopwl) RESULT(result)
    INTEGER, INTENT(IN) :: j
    REAL(4), INTENT(IN) :: wave, wl, dopwl
    REAL(4)             :: result
    REAL(4), PARAMETER :: WS(4)   = [ 6.13,  5.15,  4.24,  3.45 ]
    REAL(4), PARAMETER :: DS(4)   = [ 0.411, 0.363, 0.325, 0.293 ]
    REAL(4), PARAMETER :: ALFS(4) = [ 1.159, 1.321, 1.527, 1.783 ]
    REAL(4), PARAMETER :: FORB1(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: FORB2(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: TS(4)   = [ 5.0E3, 1.0E4, 2.0E4, 4.0E4 ]
    REAL(4), PARAMETER :: DEN = 1.0E16, DLF1 = -0.080, DLF2 = 0.0
    REAL(4) :: phihe, wtot, dtot, a_damp, sentinel

    result = 0.0
    CALL he1_isolated_profile(j, wave, wl, dopwl, WS, DS, ALFS, FORB1, FORB2, &
                             TS, DEN, DLF1, DLF2, 1.0E14, 3, sentinel)
    IF (sentinel < 0.0) THEN
      CALL read_he1_stark_tables(3, j, T(j), XNFPH(j,2), XNFHE(j,2), XNE(j), wave-wl, phihe)
      result = 1.772453 * phihe * dopwl * 10.0
      RETURN
    END IF
    CALL he1_stark_widths(j, wave, wl, WS, DS, ALFS, TS, DEN, wtot, dtot, a_damp)
    result = voigt_profile(ABS(wave - wl - dtot)/dopwl, a_damp)
  END FUNCTION he1_4387_profile

  FUNCTION he1_4921_profile(j, wave, wl, dopwl) RESULT(result)
    INTEGER, INTENT(IN) :: j
    REAL(4), INTENT(IN) :: wave, wl, dopwl
    REAL(4)             :: result
    REAL(4), PARAMETER :: WS(4)   = [ 0.002312, 0.001963, 0.001624, 0.001315 ]
    REAL(4), PARAMETER :: DS(4)   = [ 0.3932,   0.3394,   0.2950,   0.2593 ]
    REAL(4), PARAMETER :: ALFS(4) = [ 0.1207,   0.1365,   0.1564,   0.1844 ]
    REAL(4), PARAMETER :: FORB1(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: FORB2(4)= [ 0., 0., 0., 0. ]
    REAL(4), PARAMETER :: TS(4)   = [ 5.0E3, 1.0E4, 2.0E4, 4.0E4 ]
    REAL(4), PARAMETER :: DEN = 1.0E13, DLF1 = -0.150, DLF2 = 0.0
    REAL(4), PARAMETER :: FS1 = 0.0
    REAL(4) :: phihe, wtot, dtot, a_damp, wwd, sentinel

    result = 0.0
    CALL he1_isolated_profile(j, wave, wl, dopwl, WS, DS, ALFS, FORB1, FORB2, &
                             TS, DEN, DLF1, DLF2, 1.0E13, 4, sentinel)
    IF (sentinel < 0.0) THEN
      CALL read_he1_stark_tables(4, j, T(j), XNFPH(j,2), XNFHE(j,2), XNE(j), wave-wl, phihe)
      result = 1.772453 * phihe * dopwl * 10.0
      RETURN
    END IF
    CALL he1_stark_widths(j, wave, wl, WS, DS, ALFS, TS, DEN, wtot, dtot, a_damp)
    wwd    = wave - wl - dtot
    result = voigt_profile(ABS(wwd)/dopwl, a_damp)
    ! FORB1 = 0 so no forbidden component needed
  END FUNCTION he1_4921_profile


! ============================================================================
!  he1_griem_profile -- He I profile for lines in the he1_griem_profile (1974) table.
!
!  The broadening parameters for 42 He I lines are stored in the embedded
!  CHARACTER data array GRIEM0200 and parsed on first call.  The parsing
!  uses internal READ from character strings (standard F90).
!
!  For lines not in the table, falls back to a plain voigt_profile profile using
!  GAMMAR and GAMMAS.
! ============================================================================
  FUNCTION he1_griem_profile(j, wavesyn, wl, dopwl, gammar, gammas) RESULT(result)
    INTEGER, INTENT(IN) :: j
    REAL(4), INTENT(IN) :: wavesyn, wl, dopwl, gammar, gammas
    REAL(4)             :: result

    ! he1_griem_profile (1974) data: 42 lines, 5 records each (header + 4 temperatures)
    CHARACTER(LEN=55) :: GRIEM0200(210)
    DATA GRIEM0200 / &
      ' 2.00   52.2213 1s1S-4p1P 16.0                       ', &
      '   5000.  .0179     -.0112     .275     .000016      ', &
      '  10000.  .0168     -.00872    .290     .000032      ', &
      '  20000.  .0152     -.00647    .311     .000064      ', &
      '  40000.  .0135     -.00460    .341     .00013       ', &
      ' 2.00   53.7030 1s1S-3p1P 16.0                       ', &
      '   5000.  .00432    -.00277    .153     .000054      ', &
      '  10000.  .00409    -.00216    .160     .00011       ', &
      '  20000.  .00378    -.00159    .169     .00022       ', &
      '  40000.  .00341    -.00111    .183     .00043       ', &
      ' 2.00   58.4334 1s1S-2p1P 16.0                       ', &
      '   5000.  .000121   -.0000299  .012     .0024        ', &
      '  10000.  .000158   -.00000567 .010     .0049        ', &
      '  20000.  .000199    .0000220  .008     .0098        ', &
      '  40000.  .000234    .0000460  .007     .020         ', &
      ' 2.00  282.9076 2s3S-6p3P 16.0                       ', &
      '   5000.  1.79      1.07      .285      .0014        ', &
      '  10000.  1.87      .829      .276      .0027        ', &
      '  20000.  1.84      .619      .279      .0054        ', &
      '  40000.  1.72      .455      .294      .011         ', &
      ' 2.00  294.5104 2s3S-5p3P 16.0                       ', &
      '   5000.  .808      .522      .204      .0030        ', &
      '  10000.  .857      .412      .195      .0059        ', &
      '  20000.  .857      .311      .195      .012         ', &
      '  40000.  .811      .231      .203      .024         ', &
      ' 2.00  318.7746 2s3S-4p3P 16.0                       ', &
      '   5000.  .313      .219      .134      .0087        ', &
      '  10000.  .338      .176      .127      .017         ', &
      '  20000.  .344      .134      .125      .035         ', &
      '  40000.  .332      .101      .128      .069         ', &
      ' 2.00  388.8649 2s3S-3p3P 16.0                       ', &
      '   5000.  .102      .0744     .075      .050         ', &
      '  10000.  .112      .0603     .070      .099         ', &
      '  20000.  .117      .0464     .067      .20          ', &
      '  40000.  .117      .0348     .067      .40          ', &
      ' 2.00  396.4729 2s1S-4p1P 16.0                       ', &
      '   5000.  1.03      -.697      .275     .0071        ', &
      '  10000.  .996      -.504      .290     .014         ', &
      '  20000.  .877      -.374      .311     .028         ', &
      '  40000.  .776      -.266      .341     .056         ', &
      ' 2.00  402.6187 2p3P-5d3D 16.0                       ', &
      '   5000.  4.04      .541      .969      .0016        ', &
      '  10000.  3.49      .335      1.083     .0033        ', &
      '  20000.  2.96      .231      1.225     .0066        ', &
      '  40000.  2.47      .175      1.403     .013         ', &
      ' 2.00  412.0811 2p3P-5s3S 16.0                       ', &
      '   5000.  .785      .897      .171      .013         ', &
      '  10000.  .897      .894      .155      .026         ', &
      '  20000.  .984      .808      .144      .052         ', &
      '  40000.  1.01      .670      .141      .10          ', &
      ' 2.00  438.7929 2p1P-5d1D 16.0                       ', &
      '   5000.  6.13      2.52      1.159     .0017        ', &
      '  10000.  5.15      1.87      1.321     .0033        ', &
      '  20000.  4.24      1.38      1.527     .0067        ', &
      '  40000.  3.45      1.01      1.783     .013         ', &
      ' 2.00  443.7551 2p1P-5s1S 16.0                       ', &
      '   5000.  1.41      1.51       .199     .012         ', &
      '  10000.  1.57      1.43       .184     .024         ', &
      '  20000.  1.65      1.24       .177     .047         ', &
      '  40000.  1.62      .996       .179     .094         ', &
      ' 2.00  447.1488 2p3P-4d3D 16.0                       ', &
      '   5000.  1.44      .136      .589      .0058        ', &
      '  10000.  1.26      .0804     .650      .012         ', &
      '  20000.  1.09      .0614     .726      .023         ', &
      '  40000.  .927      .0549     .819      .047         ', &
      ' 2.00  471.3139 2p3P-4s3S 16.0                       ', &
      '   5000.  .342      .402      .115      .044         ', &
      '  10000.  .393      .416      .103      .088         ', &
      '  20000.  .437      .390      .095      .18          ', &
      '  40000.  .459      .335      .092      .35          ', &
      ' 2.00  492.1931 2p1P-4d1D 16.0                       ', &
      '   5000.  2.30      1.02       .683     .0061        ', &
      '  10000.  1.96      .773       .773     .012         ', &
      '  20000.  1.63      .584       .885     .024         ', &
      '  40000.  1.35      .440      1.023     .049         ', &
      ' 2.00  501.5678 2s1S-3p1P 16.0                       ', &
      '   5000.  .378      -.250      .154     .044         ', &
      '  10000.  .359      -.200      .160     .088         ', &
      '  20000.  .334      -.152      .169     .18          ', &
      '  40000.  .306      -.111      .180     .35          ', &
      ' 2.00  504.7738 2p1P-4s1S 16.0                       ', &
      '   5000.  .625       .699      .134     .038         ', &
      '  10000.  .705       .685      .123     .077         ', &
      '  20000.  .756       .611      .117     .15          ', &
      '  40000.  .760       .504      .116     .31          ', &
      ' 2.00  587.5615 2p3P-3d3D 16.0                       ', &
      '   5000.  .159      -.0881    .064      .23          ', &
      '  10000.  .170      -.0553    .061      .46          ', &
      '  20000.  .176      -.0256    .059      .92          ', &
      '  40000.  .177      -.00504   .059      1.8          ', &
      ' 2.00  667.8154 2p1P-3d1D 16.0                       ', &
      '   5000.  .423      .275       .146     .14          ', &
      '  10000.  .386      .233       .157     .27          ', &
      '  20000.  .349      .196       .169     .54          ', &
      '  40000.  .318      .161       .181     1.1          ', &
      ' 2.00  706.5176 2p3P-3s3S 16.0                       ', &
      '   5000.  .180      .215      .067      .44          ', &
      '  10000.  .207      .231      .060      .87          ', &
      '  20000.  .235      .227      .055      1.7          ', &
      '  40000.  .254      .203      .052      3.5          ', &
      ' 2.00  728.1349 2p1P-3s1S 16.0                       ', &
      '   5000.  .320       .374      .081     .33          ', &
      '  10000.  .365       .382      .073     .65          ', &
      '  20000.  .403       .355      .068     1.3          ', &
      '  40000.  .419       .303      .066     2.6          ', &
      ' 2.00  836.1694 3s3S-6p3P 16.0                       ', &
      '   5000.  15.6      9.33      .285      .035         ', &
      '  10000.  16.3      7.23      .276      .070         ', &
      '  20000.  16.1      5.38      .279      .14          ', &
      '  40000.  15.0      3.93      .293      .28          ', &
      ' 2.00  946.3596 3s3S-5p3P 16.0                       ', &
      '   5000.  8.34      5.30      .203      .099         ', &
      '  10000.  8.86      4.13      .194      .20          ', &
      '  20000.  8.88      3.06      .194      .40          ', &
      '  40000.  8.46      2.21      .201      .79          ', &
      ' 2.00  960.3418 3s1S-6p1P 16.0                       ', &
      '   5000.  40.7      -24.0     .614      .023         ', &
      '  10000.  37.3      -18.5     .656      .045         ', &
      '  20000.  33.2      -13.7     .716      .091         ', &
      '  40000.  28.8      -9.94     .796      .18          ', &
      ' 2.00 1066.7641 3p3P-6s3S 16.0                       ', &
      '   5000.  12.9      14.2      .238      .12          ', &
      '  10000.  14.7      13.9      .217      .23          ', &
      '  20000.  16.0      12.4      .203      .46          ', &
      '  40000.  16.5      10.3      .198      .92          ', &
      ' 2.00 1083.0336 2s3S-2p3P 16.0                       ', &
      '   5000.  .0493     -.0518    .028      8.3          ', &
      '  10000.  .0601     -.0557    .024      17.          ', &
      '  20000.  .0755     -.0540    .020      33.          ', &
      '  40000.  .0931     -.0465    .017      66.          ', &
      ' 2.00 1099.6561 3d3D-6p3P 16.0                       ', &
      '   5000.  27.0      16.3      .286      .08          ', &
      '  10000.  28.3      12.8      .276      .16          ', &
      '  20000.  27.9      9.59      .279      .32          ', &
      '  40000.  26.1      7.10      .294      .64          ', &
      ' 2.00 1101.3070 3s1S-5p1P 16.0                       ', &
      '   5000.  23.1      -14.3     .429      .066         ', &
      '  10000.  21.4      -11.2     .454      .13          ', &
      '  20000.  19.3      -8.49     .491      .26          ', &
      '  40000.  16.9      -6.24     .541      .53          ', &
      ' 2.00 1104.5003 3p1P-6d1D 16.0                       ', &
      '   5000.  94.5      37.7      1.742     .013         ', &
      '  10000.  79.0      27.8      1.992     .026         ', &
      '  20000.  64.8      20.4      2.310     .052         ', &
      '  40000.  52.5      14.9      2.706     .10          ', &
      ' 2.00 1196.9060 3p3P-5d3D 16.0                       ', &
      '   5000.  35.8      4.42      .967      .043         ', &
      '  10000.  31.0      2.56      1.077     .086         ', &
      '  20000.  26.4      1.65      1.215     .17          ', &
      '  40000.  22.1      1.20      1.387     .34          ', &
      ' 2.00 1252.7537 3s3S-4p3P 16.0                       ', &
      '   5000.  4.83      3.09      .130      .54          ', &
      '  10000.  5.26      2.34      .122      1.1          ', &
      '  20000.  5.45      1.65      .119      2.1          ', &
      '  40000.  5.37      1.12      .120      4.3          ', &
      ' 2.00 1284.5935 3p3P-5s3S 16.0                       ', &
      '   5000.  7.61      8.33      .165      .40          ', &
      '  10000.  8.85      8.26      .147      .81          ', &
      '  20000.  9.85      7.44      .136      1.6          ', &
      '  40000.  10.3      6.17      .132      3.2          ', &
      ' 2.00 1296.8439 3p1P-5d1D 16.0                       ', &
      '   5000.  54.3      23.1      1.149     .043         ', &
      '  10000.  45.9      17.3      1.304     .086         ', &
      '  20000.  38.1      12.8      1.501     .17          ', &
      '  40000.  31.1      9.45      1.745     .34          ', &
      ' 2.00 1298.4882 3d3D-5p3P 16.0                       ', &
      '   5000.  15.8      10.4      .205      .25          ', &
      '  10000.  16.8      8.24      .196      .51          ', &
      '  20000.  16.8      6.27      .195      1.0          ', &
      '  40000.  16.0      4.65      .203      2.0          ', &
      ' 2.00 1508.3656 3s1S-4p1P 16.0                       ', &
      '   5000.  15.1      -10.2     .277      .38          ', &
      '  10000.  14.3      -8.34     .289      .77          ', &
      '  20000.  13.2      -6.55     .307      1.5          ', &
      '  40000.  11.9      -4.95     .331      3.1          ', &
      ' 2.00 1700.2364 3p3P-4d3D 16.0                       ', &
      '   5000.  21.2      .993      .577      .32          ', &
      '  10000.  18.9      .216      .629      .64          ', &
      '  20000.  16.7      .0624     .692      1.3          ', &
      '  40000.  14.5      .129      .769      2.6          ', &
      ' 2.00 1868.5313 3d3D-4f3F 16.0                       ', &
      '   5000.  16.2      -2.77     .650      .50          ', &
      '  10000.  13.7      -1.22     .734      .99          ', &
      '  20000.  11.8      -.300     .823      2.0          ', &
      '  40000.  10.2       .131     .914      4.0          ', &
      ' 2.00 1869.7233 3d1D-4f1F 16.0                       ', &
      '   5000.  18.6      -5.39     .756      .42          ', &
      '  10000.  15.6      -3.35     .862      .84          ', &
      '  20000.  13.2      -1.95     .978      1.7          ', &
      '  40000.  11.3      -1.10     1.101     3.4          ', &
      ' 2.00 1908.9369 3p1P-4d1D 16.0                       ', &
      '   5000.  37.3      18.1      .657      .35          ', &
      '  10000.  32.3      14.0      .732      .71          ', &
      '  20000.  27.4      10.7      .828      1.4          ', &
      '  40000.  23.0      8.03      .945      2.8          ', &
      ' 2.00 1954.3172 3d3D-4p3P 16.0                       ', &
      '   5000.  12.2      8.91      .137      1.9          ', &
      '  10000.  13.2      7.25      .128      3.9          ', &
      '  20000.  13.6      5.56      .126      7.7          ', &
      '  40000.  13.3      4.11      .128      15.          ', &
      ' 2.00 2058.1299 2s1S-2p1P 16.0                       ', &
      '   5000.  .364      -.412      .038     33.          ', &
      '  10000.  .433      -.430      .033     65.          ', &
      '  20000.  .514      -.400      .029     130.         ', &
      '  40000.  .590      -.332      .026     260.         ', &
      ' 2.00 2112.0002 3p3P-4s3S 16.0                       ', &
      '   5000.  7.17      6.48      .089      4.6          ', &
      '  10000.  8.76      6.78      .077      9.2          ', &
      '  20000.  10.1      6.47      .069      18.          ', &
      '  40000.  10.9      5.62      .065      37.          ' /

    INTEGER, PARAMETER :: NLINE = 42
    INTEGER, SAVE      :: iread = 0
    REAL(4), SAVE      :: wave_t(NLINE), width_t(4,NLINE), shift_t(4,NLINE)
    REAL(4), SAVE      :: alpha_t(4,NLINE), ttab_t(4,NLINE), xnelog_t(NLINE)
    CHARACTER(LEN=10), SAVE :: trans_t(NLINE)

    INTEGER :: il, i, it
    REAL(4) :: e, temp, x, xx, w_damp, d_ratio, alf, xnfhp, xnfhep
    REAL(4) :: vm1, rhom, sigma, wtot, dtot, a_damp, v
    REAL(4) :: dummy_code, dummy_xnelog
    CHARACTER(LEN=55) :: cbuf  ! scratch for internal READ from PARAMETER array

    ! --- Parse table on first call ---
    IF (iread == 0) THEN
      DO il = 1, NLINE
        cbuf = GRIEM0200(il*5-4)
        READ(cbuf, '(F5.2,F10.4,A10,F5.1)') &
             dummy_code, wave_t(il), trans_t(il), xnelog_t(il)
        DO i = 1, 4
          cbuf = GRIEM0200(il*5-4+i)
          READ(cbuf, '(F10.0,4F10.6)') &
               ttab_t(i,il), width_t(i,il), shift_t(i,il), alpha_t(i,il), dummy_xnelog
          shift_t(i,il) = shift_t(i,il) / width_t(i,il)
        END DO
      END DO
      iread = 1
    END IF

    ! --- Search for matching line ---
    il = 0
    DO i = 1, NLINE
      IF (ABS(wl - wave_t(i)) < 0.1) THEN
        il = i
        EXIT
      END IF
    END DO

    ! Not in table: plain voigt_profile with GAMMAR + GAMMAS
    IF (il == 0) THEN
      v      = ABS(wavesyn - wl) / dopwl
      a_damp = (gammar + gammas*XNE(j)) / (dopwl/wl)
      result = voigt_profile(v, a_damp)
      RETURN
    END IF

    ! In table: he1_griem_profile Stark broadening
    e      = XNE(j)
    temp   = MAX(MIN(T(j), 8.0E4), 5.0E3)
    xnfhp  = XNFPH(j,2)
    xnfhep = XNFHE(j,2)

    it = 2
    DO WHILE (it < 4 .AND. ttab_t(it,il) <= temp)
      it = it + 1
    END DO
    x  = (temp - ttab_t(it-1,il)) / (ttab_t(it,il) - ttab_t(it-1,il))
    xx = e / (10.0**xnelog_t(il))

    w_damp  = xx * (x*width_t(it,il) + (1.0-x)*width_t(it-1,il))
    d_ratio = x*shift_t(it,il) + (1.0-x)*shift_t(it-1,il)
    alf     = xx**0.25 * (x*alpha_t(it,il) + (1.0-x)*alpha_t(it-1,il))
    xx      = xnfhp / e
    vm1     = 8.78 * (xx + 2.0*(1.0-xx)) / SQRT(temp)
    rhom    = 1.0 / (4.19*e)**(1.0/3.0)
    sigma   = 1.885E14 * w_damp * rhom * vm1 / (wl*10.0)**2
    x       = alf**(8.0/9.0) / sigma**(1.0/3.0)
    wtot    = w_damp * (1.0 + 1.36*x) * 0.1
    dtot    = w_damp * d_ratio * (1.0 + 2.36*x/ABS(d_ratio)) * 0.1
    ! gammar damping: use dopwl/wl as Doppler unit for frequency-based ratio
    a_damp  = wtot/dopwl + gammar/(dopwl/wl)
    result  = voigt_profile(ABS(wavesyn - wl - dtot)/dopwl, a_damp)

  END FUNCTION he1_griem_profile


! ============================================================================
!  he1_dimitri_profile -- He I profile with he1_dimitri_profile (Mihalas et al.) broadening tables.
!
!  Similar to he1_griem_profile but uses a 7-line embedded table with separate electron,
!  proton, and He broadening widths and shifts.  Parsed on first call.
! ============================================================================
  FUNCTION he1_dimitri_profile(j, wavesyn, wl, dopwl, gammar, gammas) RESULT(result)
    INTEGER, INTENT(IN) :: j
    REAL(4), INTENT(IN) :: wavesyn, wl, dopwl, gammar, gammas
    REAL(4)             :: result

    CHARACTER(LEN=61) :: HE1DATTAB(35)
    DATA HE1DATTAB / &
      ' 2.00  381.9624 2p3P-6d3D 13.0                               ', &
      '  5000. 1.47E-02-2.04E-04 9.31E-03 7.87E-03 0.00E+00 0.00E+00', &
      ' 10000. 1.29E-02-4.41E-05 1.06E-02 9.31E-03 8.90E-03 7.72E-03', &
      ' 20000. 1.10E-02-1.37E-04 1.72E-02 1.13E-02 1.02E-02 8.99E-03', &
      ' 40000. 9.33E-03-2.31E-05 1.08E-02 1.32E-02 1.13E-02 1.06E-02', &
      ' 2.00  386.7494 2p3P-6s3S 16.0                               ', &
      '  5000. 2.74E+00 1.79E+00 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 10000. 2.91E+00 1.86E+00 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 20000. 2.85E+00 1.51E+00 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 40000. 3.10E+00 1.10E+00 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 2.00  392.6544 2p1P-6d1D 13.0                               ', &
      '  5000  7.07E-02 7.13E-03 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 10000. 6.07E-02 4.27E-03 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 20000. 5.09E-02 1.96E-03 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 40000. 4.19E-02 5.60E-04 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 2.00  400.9256 2p1P-7d1D 13.0                               ', &
      '  5000  3.96E-02 6.23E-03 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 10000. 3.42E-02 2.78E-03 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 20000. 2.88E-02 1.85E-03 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 40000. 2.38E-02 1.07E-03 3.05E-02 3.67E-02 0.00E+00 0.00E+00', &
      ' 2.00  402.40 2p1P-7s1S 13.0                                 ', &
      '  5000  9.24E-03 6.25E-03 2.09E-03 1.91E-03 1.79E-03 1.63E-03', &
      ' 10000. 8.80E-03 5.30E-03 2.34E-03 2.17E-03 2.01E-03 1.86E-03', &
      ' 20000. 9.18E-03 3.77E-03 2.63E-03 2.46E-03 2.42E-03 2.72E-03', &
      ' 40000. 9.19E-03 2.32E-03 2.96E-03 2.77E-03 2.54E-03 2.38E-03', &
      ' 2.00  414.3759 2p1P-6d1D 13.0                               ', &
      '  5000. 2.21E-02 3.63E-03 0.00E+00 0.00E+00 0.00E+00 0.00E+00', &
      ' 10000. 1.90E-02 2.37E-03 1.85E-02 1.67E-02 0.00E+00 0.00E+00', &
      ' 20000. 1.59E-02 1.12E-03 1.87E-02 1.96E-02 0.00E+00 0.00E+00', &
      ' 40000. 1.31E-02 6.75E-04 1.65E-02 2.16E-02 1.89E-02 1.90E-02', &
      ' 2.00  416.8972 2p1P-6s1S 13.0                               ', &
      '  5000. 4.83E-03 3.41E-03 1.08E-03 1.00E-03 9.29E-04 8.57E-04', &
      ' 10000. 4.72E-03 3.06E-03 1.21E-03 1.13E-03 1.04E-03 9.72E-04', &
      ' 20000. 4.84E-03 2.27E-03 1.36E-03 1.28E-03 1.17E-03 1.10E-03', &
      ' 40000. 4.97E-03 1.55E-03 1.53E-03 1.44E-03 1.32E-03 1.24E-03' /

    INTEGER, PARAMETER :: NLINE = 7
    INTEGER, SAVE      :: iread = 0
    REAL(4), SAVE      :: wave_t(NLINE), xnelog_t(NLINE), ttab_t(4,NLINE)
    REAL(4), SAVE      :: width_t(4,NLINE),  shift_t(4,NLINE)
    REAL(4), SAVE      :: widthp_t(4,NLINE), shiftp_t(4,NLINE)
    REAL(4), SAVE      :: widthhe_t(4,NLINE),shifthe_t(4,NLINE)
    CHARACTER(LEN=10), SAVE :: trans_t(NLINE)

    INTEGER :: il, i, it
    REAL(4) :: e, temp, x, xx, xxh, xxhe, w_e, w_p, w_he
    REAL(4) :: d_e, d_p, d_he, wtot, dtot, a_damp, v
    REAL(4) :: dummy_code
    REAL(4) :: xnfhp, xnfhep
    CHARACTER(LEN=61) :: cbuf  ! scratch for internal READ from PARAMETER array

    ! Parse table on first call
    IF (iread == 0) THEN
      DO il = 1, NLINE
        cbuf = HE1DATTAB(il*5-4)
        READ(cbuf, '(F5.2,F10.4,A10,F5.1)') &
             dummy_code, wave_t(il), trans_t(il), xnelog_t(il)
        DO i = 1, 4
          cbuf = HE1DATTAB(il*5-4+i)
          READ(cbuf, '(F7.0,6E9.2)') &
               ttab_t(i,il), width_t(i,il), shift_t(i,il), &
               widthp_t(i,il), shiftp_t(i,il), widthhe_t(i,il), shifthe_t(i,il)
          shift_t(i,il)  = shift_t(i,il)  / width_t(i,il)
          IF (widthp_t(i,il)  > 0.0) shiftp_t(i,il)  = shiftp_t(i,il)  / widthp_t(i,il)
          IF (widthhe_t(i,il) > 0.0) shifthe_t(i,il) = shifthe_t(i,il) / widthhe_t(i,il)
        END DO
      END DO
      iread = 1
    END IF

    ! Search for matching line
    il = 0
    DO i = 1, NLINE
      IF (ABS(wl - wave_t(i)) < 0.1) THEN
        il = i
        EXIT
      END IF
    END DO

    IF (il == 0) THEN
      v      = ABS(wavesyn - wl) / dopwl
      a_damp = (gammar + gammas*XNE(j)) / (dopwl/wl)
      result = voigt_profile(v, a_damp)
      RETURN
    END IF

    e      = XNE(j)
    temp   = MAX(MIN(T(j), 8.0E4), 5.0E3)
    xnfhp  = XNFPH(j,2)
    xnfhep = XNFHE(j,2)

    it = 2
    DO WHILE (it < 4 .AND. ttab_t(it,il) <= temp)
      it = it + 1
    END DO
    x   = (temp - ttab_t(it-1,il)) / (ttab_t(it,il) - ttab_t(it-1,il))
    xx  = e    / (10.0**xnelog_t(il))
    xxh = xnfhp  / (10.0**xnelog_t(il))
    xxhe= xnfhep / (10.0**xnelog_t(il))

    w_e  = xx   * (x*width_t(it,il)   + (1.0-x)*width_t(it-1,il))
    w_p  = xxh  * (x*widthp_t(it,il)  + (1.0-x)*widthp_t(it-1,il))
    w_he = xxhe * (x*widthhe_t(it,il) + (1.0-x)*widthhe_t(it-1,il))
    d_e  = shift_t(it,il)  + (1.0-x)*shift_t(it-1,il)    ! note: original typo preserved
    d_p  = x*shiftp_t(it,il)  + (1.0-x)*shiftp_t(it-1,il)
    d_he = x*shifthe_t(it,il) + (1.0-x)*shifthe_t(it-1,il)

    wtot   = (w_e + w_p + w_he) * 0.1 / 2.0   ! Angstrom -> nm, then half-half-width
    dtot   = (w_e*d_e + w_p*d_p + w_he*d_he) * 0.1
    a_damp = wtot/dopwl + gammar/(dopwl/wl)
    result = voigt_profile(ABS(wavesyn - wl - dtot)/dopwl, a_damp)

  END FUNCTION he1_dimitri_profile


! ============================================================================
!  read_he1_stark_tables -- read and interpolate BCS (Barnard, Cooper, Saraph) tabulated
!             He I profiles for lines 4471, 4026, 4387, 4921 nm.
!
!  Reads binary table data from unit 18 on first call (initialised by
!  the sentinel value DLAM(1,1) == -150.0).  Subsequent calls interpolate
!  in temperature, electron density, and wavelength offset.
!
!  Arguments:
!    LINE   -- line index: 1=4471, 2=4026, 3=4387, 4=4921
!    J      -- atmospheric depth index
!    TEMP   -- temperature (K)
!    XNFHP  -- proton number density (cm^{-3})
!    XNFHEP -- He+ number density (cm^{-3})
!    XNE_IN -- electron number density (cm^{-3})
!    DLNM   -- wavelength offset from line centre (nm)
!    PHIHE  -- output: normalised line profile value
! ============================================================================
  SUBROUTINE read_he1_stark_tables(line, j, temp, xnfhp, xnfhep, xne_in, dlnm, phihe)
    INTEGER, INTENT(IN)  :: line, j
    REAL(4), INTENT(IN)  :: temp, xnfhp, xnfhep, xne_in, dlnm
    REAL(4), INTENT(OUT) :: phihe

    INTEGER, PARAMETER :: NDLAM(4) = [ 142, 196, 204, 142 ]
    INTEGER, PARAMETER :: NXNE(4)  = [ 7,     8,   8,   7 ]
    REAL(4), PARAMETER :: XNE1(4)  = [ 13.0, 14.0, 14.0, 13.0 ]

    ! Tables read from fort.18; sentinel DLAM(1,1) = 0 means not yet read
    REAL(4), SAVE :: dlam(204,4)
    REAL(4), SAVE :: phihp(4,7,142), phihep(4,7,142)
    REAL(4), SAVE :: phi4026(4,8,196), phi4387(4,8,204)
    REAL(4), SAVE :: phi4921(4,7,142)
    INTEGER, SAVE :: itemp1_bcs = 0
    INTEGER, SAVE :: jsave = 0
    REAL(4), SAVE :: philam(204)
    LOGICAL, SAVE :: initialised = .FALSE.

    CHARACTER(LEN=8) :: title1, title2
    INTEGER :: il, ne, it, ip, i
    REAL(4) :: at, bt, wt_t, ap, bp, wt_p, c1w1w, c1ww, cw1w, cww
    REAL(4) :: fne, dwl, phi_tmp(8), xxh, xxhe, phinorm, dl, a, b
    INTEGER :: it_loc, ip_loc

    ! --- Read table on first call ---
    IF (.NOT. initialised) THEN
      dlam = 0.0

      OPEN(UNIT=18, FILE=trim(DATADIR)//'he1tables.dat', FORM='FORMATTED', STATUS='OLD')

      ! 4471: 142 wavelength points, 7 Ne values, 8 entries (4 T x 2 species)
      READ(18, '(A8)') title1;  READ(18, '(A8)') title2
      DO il = 1, 142
        DO ne = 1, 7
          READ(18, '(1X,F5.1,F8.2,8F7.3)') fne, dwl, (phi_tmp(i), i=1,8)
          DO it = 1, 4
            phihp(it,ne,il)  = phi_tmp(it)
            phihep(it,ne,il) = phi_tmp(it+4)
          END DO
        END DO
        dlam(il,1) = dwl - 150.0
      END DO

      ! 4026: 196 wavelength points, 8 Ne values, 4 T entries
      READ(18, '(A8)') title1;  READ(18, '(A8)') title2
      DO il = 1, 196
        DO ne = 1, 8
          READ(18, '(1X,F5.1,F8.2,4F7.3)') fne, dwl, (phi4026(it,ne,il), it=1,4)
        END DO
        dlam(il,2) = dwl - 150.0
      END DO

      ! 4387: 204 wavelength points, 8 Ne values, 4 T entries
      READ(18, '(A8)') title1;  READ(18, '(A8)') title2
      DO il = 1, 204
        DO ne = 1, 8
          READ(18, '(1X,F5.1,F8.2,8F7.3)') fne, dwl, (phi_tmp(i), i=1,8)
          DO it = 1, 4
            phi4387(it,ne,il) = phi_tmp(it)
          END DO
        END DO
        dlam(il,3) = dwl - 150.0
      END DO

      ! 4921: 142 wavelength points, 7 Ne values, 8 entries (4 T x 2 species)
      READ(18, '(A8)') title1;  READ(18, '(A8)') title2
      DO il = 1, 142
        DO ne = 1, 7
          READ(18, '(1X,F5.1,F8.2,8F7.3)') fne, dwl, (phi_tmp(i), i=1,8)
          DO it = 1, 4
            phihp(it,ne,il)  = phi_tmp(it)
            phihep(it,ne,il) = phi_tmp(it+4)
          END DO
        END DO
        dlam(il,4) = dwl - 150.0
      END DO

      CLOSE(UNIT=18)
      initialised = .TRUE.
    END IF

    ! --- Interpolate in T and Ne if depth point changed ---
    IF (j * line /= jsave) THEN
      ! Temperature interpolation (4 points: log T = 3.699, 4.000, 4.301, 4.602)
      at   = LOG10(temp)
      bt   = (at - 3.698970) / 0.3010300 + 1.0
      it_loc = MAX(MIN(INT(bt + 1.0E-5), 3), 1)
      wt_t = bt - it_loc

      ! Electron density interpolation
      ap   = LOG10(xne_in)
      ap   = MAX(XNE1(line), ap)
      bp   = (ap - XNE1(line)) / 0.5 + 1.0
      ip_loc = MAX(MIN(INT(bp + 1.0E-5), NXNE(line)-1), 1)
      wt_p = bp - ip_loc

      c1w1w = (1.0-wt_p)*(1.0-wt_t)
      c1ww  = (1.0-wt_p)*wt_t
      cw1w  = wt_p*(1.0-wt_t)
      cww   = wt_p*wt_t

      xxh  = xnfhp  / xne_in
      xxhe = xnfhep / xne_in

      SELECT CASE (line)
        CASE (1)   ! 4471
          DO i = 1, NDLAM(line)
            philam(i) = xxh * 10.0**( c1w1w*phihp(it_loc  ,ip_loc  ,i) + &
                                       c1ww *phihp(it_loc+1,ip_loc  ,i) + &
                                       cw1w *phihp(it_loc  ,ip_loc+1,i) + &
                                       cww  *phihp(it_loc+1,ip_loc+1,i) ) + &
                         xxhe* 10.0**( c1w1w*phihep(it_loc  ,ip_loc  ,i) + &
                                       c1ww *phihep(it_loc+1,ip_loc  ,i) + &
                                       cw1w *phihep(it_loc  ,ip_loc+1,i) + &
                                       cww  *phihep(it_loc+1,ip_loc+1,i) )
          END DO
        CASE (2)   ! 4026
          DO i = 1, NDLAM(line)
            philam(i) = 10.0**( c1w1w*phi4026(it_loc  ,ip_loc  ,i) + &
                                 c1ww *phi4026(it_loc+1,ip_loc  ,i) + &
                                 cw1w *phi4026(it_loc  ,ip_loc+1,i) + &
                                 cww  *phi4026(it_loc+1,ip_loc+1,i) )
          END DO
        CASE (3)   ! 4387
          DO i = 1, NDLAM(line)
            philam(i) = 10.0**( c1w1w*phi4387(it_loc  ,ip_loc  ,i) + &
                                 c1ww *phi4387(it_loc+1,ip_loc  ,i) + &
                                 cw1w *phi4387(it_loc  ,ip_loc+1,i) + &
                                 cww  *phi4387(it_loc+1,ip_loc+1,i) )
          END DO
        CASE (4)   ! 4921 (same table structure as 4471)
          DO i = 1, NDLAM(line)
            philam(i) = xxh * 10.0**( c1w1w*phihp(it_loc  ,ip_loc  ,i) + &
                                       c1ww *phihp(it_loc+1,ip_loc  ,i) + &
                                       cw1w *phihp(it_loc  ,ip_loc+1,i) + &
                                       cww  *phihp(it_loc+1,ip_loc+1,i) ) + &
                         xxhe* 10.0**( c1w1w*phihep(it_loc  ,ip_loc  ,i) + &
                                       c1ww *phihep(it_loc+1,ip_loc  ,i) + &
                                       cw1w *phihep(it_loc  ,ip_loc+1,i) + &
                                       cww  *phihep(it_loc+1,ip_loc+1,i) )
          END DO
      END SELECT

      ! Normalise
      CALL trapz_integrate(dlam(1:NDLAM(line),line), philam(1:NDLAM(line)), phinorm, NDLAM(line), 0.0)
      DO i = 1, NDLAM(line)
        philam(i) = LOG10(philam(i) / phinorm)
      END DO
      jsave = j * line
    END IF

    ! --- Interpolate in wavelength offset ---
    dl = dlnm * 10.0   ! nm -> Angstrom
    phihe = 0.0
    DO i = 2, NDLAM(line)
      IF (dl <= dlam(i,line) .OR. i == NDLAM(line)) THEN
        a     = (dlam(i,line) - dl) / (dlam(i,line) - dlam(i-1,line))
        b     = (dl - dlam(i-1,line)) / (dlam(i,line) - dlam(i-1,line))
        phihe = 10.0**(a*philam(i-1) + b*philam(i))
        RETURN
      END IF
    END DO

  END SUBROUTINE read_he1_stark_tables


! ============================================================================
!  compute_line_opacity -- add NLTE line opacities (from TAPE19/unit 19) to BUFFER for
!            depth point J.
!
!  This handles five line types dispatched on the integer TYPE field:
!    TYPE =  0   normal metal line (voigt_profile profile)
!    TYPE =  3   PRD line          (treated as TYPE=0)
!    TYPE = -1   hydrogen line (H I)
!    TYPE = -2   deuterium line    (treated as H I with heavier Doppler width)
!    TYPE =  1   autoionizing (Shore) line
!    TYPE =  2   coronal line      (skipped)
!    TYPE < -2   He I line
!  Any other TYPE >= 4 is treated as a merged-continuum edge.
!
!  F77 EQUIVALENCE aliases resolved as local variables:
!    CGF   = GF   (oscillator strength, also G)
!    ASHORE= GAMMAS  (asymmetry parameter q for autoionizing profile)
!    BSHORE= GAMMAW  (width parameter for autoionizing profile)
!    XSECT = GAMMAR  (also GAUNT)
!    NLAST = TYPE    (principal quantum number of series limit)
!
!  Arguments:
!    J        -- current depth index (1..NRHOX)
!    N19      -- number of NLTE lines on unit 19
!    CUTOFF   -- minimum kappa/continuum ratio to include a line
!    VELSHIFT -- radial velocity shift for this depth (km/s), REAL(4)
!    IFVAC    -- 1 = wavelengths already in vacuum, 0 = convert air->vacuum
!    LINOUT   -- >= 0: write (ILINE, KAPCEN) records to unit 15
! ============================================================================
  SUBROUTINE compute_line_opacity(j, n19, cutoff, velshift, ifvac, linout)

    INTEGER,  INTENT(IN) :: j, n19, ifvac, linout
    REAL(4),  INTENT(IN) :: cutoff, velshift

    ! ----- local scalars -----
    REAL(4)  :: bolt, bolth, oldelo, oldeloh
    REAL(4)  :: kappa0, kappa, kapcen, kapmin
    REAL(4)  :: kappa0red, kappared, kappa0blue, kappablue
    REAL(4)  :: adamp, dopwl, vvoigt, epsil, frelin, freq
    REAL(4)  :: xsectg, tail, dnbuff
    REAL(4)  :: edgeblue
    REAL(4)  :: v2, nelem_r
    REAL(4)  :: hfac(kw), hefac(kw), h2fac(kw)
    REAL(4)  :: alpha_bao            ! Barklem-Anstee-O'Mara exponent
    INTEGER  :: iline, ncon, nelionx, itype, nlast_loc
    INTEGER  :: nbuff, nbuff1, nbuff2, nbuff3, ibuff
    INTEGER  :: minred, maxblue, ixwl, nelem_i
    INTEGER  :: i, k, n
    REAL(8)  :: wave, wcon, wmerge, wshift, wtail
    REAL(8)  :: bluecut, wlplus1, wlplus2, redcut, wlminus1, wlminus2
    REAL(8)  :: dopratio, wl_loc
    REAL(8)  :: emergeh_loc(kw)
    REAL(4)  :: dopph(kw)
    ! local aliases replacing F77 EQUIVALENCE on GAMMAR/GAMMAS/GAMMAW/GF
    REAL(4)  :: cgf_loc, ashore_loc, bshore_loc, xsect_loc
    REAL(4)  :: gammar_loc, gammas_loc, gammaw_loc
    REAL(4)  :: elo_loc, gf_loc
    REAL(4)  :: nstark, ndopp, nmerge, inglis  ! declared REAL in F77

    ! ----- saved between calls -----
    INTEGER, SAVE :: itemp1 = 0
    REAL(8), SAVE :: ehyd(100)    = 0.0D0
    REAL(8), SAVE :: alphahyd(99) = 0.0D0

    ! ----- DATA: hydrogen series limits (cm^{-1}) used for merged continuum -----
    REAL(8), PARAMETER :: conth(15) = [ &
      109678.764D0, 27419.659D0, 12186.462D0,  6854.871D0, 4387.113D0, &
        3046.604D0,  2238.320D0,  1713.711D0,  1354.044D0, 1096.776D0, &
         906.426D0,   761.650D0,   648.980D0,   559.579D0,  487.456D0 ]

    ! ----- DATA: CONTX(26,17) -- ionisation edge wavenumbers -----
    ! Column-major order from the F77 DATA statement (26 rows, 17 columns).
    ! N*val repeat syntax is invalid in F90 array constructors; zeros spelled out.
    REAL(8), PARAMETER :: contx(26,17) = RESHAPE( [ &
      ! col 1  -- H I (1.00)
      109678.764D0, 27419.659D0, 12186.462D0,  6854.871D0, 4387.113D0, &
        3046.604D0,  2238.320D0,  1713.711D0,  1354.044D0, 1096.776D0, &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,            &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! col 2  -- He I (2.00)
      198310.760D0, 38454.691D0, 32033.214D0, 29223.753D0, 27175.760D0, &
       15073.868D0, 0.D0,0.D0,0.D0,0.D0,                                &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! col 3  -- He II (2.01)
      438908.850D0,109726.529D0, 48766.491D0, 27430.925D0, 17555.715D0, &
       12191.437D0, 0.D0,0.D0,0.D0,0.D0,                                &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! col 4  -- C I (6.00)
       90883.840D0, 90867.420D0, 90840.420D0, 90820.420D0, 90804.000D0, &
       90777.000D0, 80691.180D0, 80627.760D0, 69235.820D0, 69172.400D0, &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! col 5  -- C II (6.01)
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! col 6  -- Mg I (12.00)
       61671.020D0, 39820.615D0, 39800.556D0, 39759.842D0,               &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,                                                        &
      ! col 7  -- Mg II (12.01)
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! col 8  -- Al I (13.00)
       48278.370D0, 48166.309D0,                                          &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,                                             &
      ! col 9  -- Al II (13.01)
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! col 10 -- Si I (14.00)
       66035.000D0, 65957.885D0, 65811.843D0, 65747.550D0, 65670.435D0, &
       65524.393D0, 59736.150D0, 59448.700D0, 50640.630D0, 50553.180D0, &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
      ! cols 11-17 -- all zeros (Si II, Ca I, Ca II, O I, O II, Na I, K I)
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,                                  &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,             &
        0.D0,0.D0,0.D0,0.D0,0.D0,0.D0                                   &
      ], SHAPE=[26,17] )

    ! ------------------------------------------------------------------
    !  One-time initialisation: EHYD, ALPHAHYD, EMERGE, EMERGEH
    !  Triggered whenever ITEMP changes (new temperature structure).
    ! ------------------------------------------------------------------
    IF (itemp /= itemp1) THEN
      ehyd(1) = 0.0D0
      ehyd(2) = 82259.105D0
      ehyd(3) = 97492.302D0
      ehyd(4) = 102823.893D0
      ehyd(5) = 105291.651D0
      ehyd(6) = 106632.160D0
      ehyd(7) = 107440.444D0
      ehyd(8) = 107965.051D0
      DO n = 9, 100
        ehyd(n) = 109678.764D0 - 109677.576D0 / DBLE(n)**2
      END DO
      DO n = 1, 99
        alphahyd(n) = 1.0D7 / (ehyd(n+1) - ehyd(n))
      END DO
      DO k = 1, nrhox
        ! Inglis-Teller merger level (empirical)
        inglis = 1600.0 / xne(k)**(2.0/15.0)
        nmerge = inglis - 1.5
        emerge(k)     = 109737.312D0 / DBLE(nmerge)**2
        emergeh_loc(k)= 109677.576D0 / DBLE(nmerge)**2
      END DO
      itemp1 = itemp
    END IF

    ! ------------------------------------------------------------------
    !  Per-call initialisation
    ! ------------------------------------------------------------------
    bolt      = 1.0
    bolth     = 1.0
    oldelo    = 1.0E30 * REAL(j * itemp)
    oldeloh   = 1.0E30 * REAL(j * itemp)
    dopratio  = 1.0D0 + DBLE(velshift) / CLIGHT_KMS
    alpha_bao = 0.0

    REWIND 19
    line_loop: DO iline = 1, n19

      READ(19) wl_loc, elo_loc, gf_loc, nblo, nbup, nelion, itype, ncon, nelionx, &
               gammar_loc, gammas_loc, gammaw_loc, nbuff
      ! (Barklem variant would also read alpha_bao before nbuff)

      ! Resolve F77 EQUIVALENCE aliases
      cgf_loc   = gf_loc        ! GF = G = CGF
      ashore_loc= gammas_loc    ! GAMMAS = ASHORE
      bshore_loc= gammaw_loc    ! GAMMAW = BSHORE
      xsect_loc = gammar_loc    ! GAMMAR = XSECT = GAUNT
      nlast_loc = itype         ! TYPE = NLAST (for merged-continuum branch)

      wl_loc = wl_loc * dopratio

      ! ------ Dispatch on line type ------
      SELECT CASE (itype)
        CASE (2)
          CYCLE line_loop          ! coronal -- skip
        CASE (0)
          GOTO 200                 ! normal metal line
        CASE (-1, -2)
          GOTO 600                 ! hydrogen / deuterium
        CASE (1)
          GOTO 700                 ! autoionizing (Shore)
        CASE (3)
          GOTO 200                 ! PRD -- treat as normal
        CASE (:-3)
          GOTO 200                 ! He I not yet implemented -- treat as normal
        CASE DEFAULT
          ! TYPE >= 4: merged continuum edge
      END SELECT

      ! ------------------------------------------------------------------
      !  MERGED CONTINUUM EDGE
      ! ------------------------------------------------------------------
      wshift = 1.0D7 / (1.0D7/wl_loc - 109737.312D0/DBLE(nlast_loc)**2)
      wmerge = 1.0D7 / (1.0D7/wl_loc - emerge(j))
      IF (nelion == 1) THEN
        wshift = 1.0D7 / (1.0D7/wl_loc - 109677.576D0/DBLE(nlast_loc)**2)
        wmerge = 1.0D7 / (1.0D7/wl_loc - emergeh_loc(j))
      END IF
      IF (wmerge < 0.0D0) wmerge = wshift + wshift
      wmerge = MAX(wmerge, wshift)
      wmerge = MIN(wshift + wshift, wmerge)
      wtail  = 1.0D7 / (1.0D7/wmerge - 500.0D0)
      IF (wtail < 0.0D0) wtail = wmerge + wmerge
      wtail  = MIN(wmerge + wmerge, wtail)
      IF (ifvac == 0) THEN
        wmerge = vac_to_air(wmerge) * dopratio
        wtail  = vac_to_air(wtail)  * dopratio
      END IF
      ixwl   = INT(LOG(wl_loc)  / ratiolg)
      edgeblue = EXP(REAL(ixwl) * REAL(ratiolg))
      IF (edgeblue > REAL(wl_loc)) ixwl = ixwl - 1
      nbuff1 = ixwl + 1 - ixwlbeg + 1
      ixwl   = INT(LOG(wmerge) / ratiolg + 0.5D0)
      nbuff2 = ixwl - ixwlbeg + 1
      ixwl   = INT(LOG(wtail)  / ratiolg + 0.5D0)
      nbuff3 = ixwl - ixwlbeg + 1
      IF (nbuff1 > length) CYCLE line_loop
      IF (nbuff3 < 1)      CYCLE line_loop
      dnbuff = REAL(nbuff3 - nbuff2)
      nbuff1 = MAX(nbuff1, 1)
      xsectg = gf_loc
      kappa  = xsectg * REAL(xnfpel(nelion)) * REAL(EXP(-elo_loc * hckt(j)))
      tail   = 1.0
      DO ibuff = nbuff1, MIN(nbuff3, length)
        IF (ibuff > nbuff2) tail = REAL(nbuff3 - ibuff) / dnbuff
        buffer(ibuff) = buffer(ibuff) + kappa * tail
      END DO
      IF (linout >= 0) WRITE(15) iline, kappa
      mlines = mlines + 1
      CYCLE line_loop

      ! ------------------------------------------------------------------
200   CONTINUE   ! NORMAL METAL LINE (TYPE 0, 3, or He I fallthrough)
      ! ------------------------------------------------------------------
      kappa0 = REAL(cgf_loc * xnfdop(nelion))   ! xnfdop now R8; result narrowed to R4 kappa0
      kapmin = continuum(MIN(MAX(nbuff,1),length)) * cutoff
      IF (kappa0 < kapmin) CYCLE line_loop
      IF (elo_loc /= oldelo) THEN
        bolt   = REAL(EXP(-elo_loc * hckt(j)))   ! Boltzmann factor in R8, narrowed to R4
        oldelo = elo_loc
      END IF
      kappa0 = kappa0 * bolt
      IF (kappa0 < kapmin) CYCLE line_loop
      mlines = mlines + 1
      wcon  = 0.0D0
      wtail = 0.0D0
      IF (ncon > 0) THEN
        wcon  = 1.0D7 / (contx(ncon,nelionx) - emerge(j)) * dopratio
        wtail = 1.0D7 / (contx(ncon,nelionx) - emerge(j) - 500.0D0) * dopratio
      END IF
      ! Barklem-Anstee-O'Mara van der Waals (when alpha_bao /= 0)
      IF (alpha_bao /= 0.0) THEN
        nelem_i = INT(nelion/6) + 1
        v2 = (1.0 - alpha_bao) / 2.0
        hfac(j)  = (t(j)/10000.0)**v2
        hefac(j) = 0.628 * (2.0991E-4*t(j)*(0.25 + 1.008/atmass(nelem_i)))**v2
        h2fac(j) = 1.08  * (2.0991E-4*t(j)*(0.50 + 1.008/atmass(nelem_i)))**v2
        txnxn(j) = xnfh(j)*hfac(j) + xnfhe(j,1)*hefac(j) + xnfh2(j)*h2fac(j)
      END IF
      adamp  = REAL((gammar_loc + gammas_loc*xne(j) + gammaw_loc*txnxn(j)) / dopple(nelion))
      kapcen = kappa0 * voigt_profile(0.0, adamp)
      IF (linout >= 0) WRITE(15) iline, kapcen
      dopwl  = REAL(dopple(nelion)) * REAL(wl_loc)
      ! Red wing
      IF (wl_loc <= wlend) THEN
        minred = MAX(1, nbuff)
        wave   = wbegin * ratio**(minred-1)
        DO ibuff = minred, length
          IF (wave >= wcon) THEN
            vvoigt = ABS(REAL(wave - wl_loc)) / dopwl
            kappa  = kappa0 * voigt_profile(vvoigt, adamp)
            IF (wave < wtail) kappa = kappa * REAL((wave-wcon)/(wtail-wcon))
            IF (kappa < continuum(ibuff)*cutoff) EXIT
            buffer(ibuff) = buffer(ibuff) + kappa
          END IF
          wave = wave * ratio
        END DO
        IF (minred == 1) CYCLE line_loop
        IF (wl_loc < wbegin) CYCLE line_loop
      END IF
      ! Blue wing
      ibuff    = MIN(length+1, nbuff)
      maxblue  = ibuff - 1
      wave     = wbegin * ratio**(ibuff-1)
      DO i = 1, maxblue
        ibuff  = ibuff - 1
        wave   = wave / ratio
        IF (wave < wcon) CYCLE
        vvoigt = ABS(REAL(wave - wl_loc)) / dopwl
        kappa  = kappa0 * voigt_profile(vvoigt, adamp)
        IF (wave < wtail) kappa = kappa * REAL((wave-wcon)/(wtail-wcon))
        buffer(ibuff) = buffer(ibuff) + kappa
        IF (kappa < continuum(ibuff)*cutoff) CYCLE line_loop
      END DO
      CYCLE line_loop

      ! ------------------------------------------------------------------
600   CONTINUE   ! HYDROGEN / DEUTERIUM LINE
      ! ------------------------------------------------------------------
      kappa0 = REAL(cgf_loc * xnfdop(1))
      kapmin = continuum(MIN(MAX(nbuff,1),length)) * cutoff
      IF (nbup == 2) kapmin = continuum(length) * cutoff
      IF (kappa0 < kapmin) CYCLE line_loop
      bolth    = REAL(EXP(-elo_loc * hckt(j)))
      oldeloh  = elo_loc
      dopph(j) = REAL(dopple(1))
      IF (itype == -2) dopph(j) = dopph(j) / 1.4142    ! deuterium
      kappa0 = kappa0 * bolth
      IF (kappa0 < kapmin) CYCLE line_loop
      IF (linout >= 0) WRITE(15) iline, kappa0
      mlines = mlines + 1
      ! Alpha (Nbup=Nblo+1) and beta-blue (Nbup=Nblo+2) treated as isolated
      IF (ncon == 0)            GOTO 620
      IF (nbup == nblo+1)       GOTO 620
      IF (nbup == nblo+2)       GOTO 630
      ! General Balmer/Paschen/... with merged continuum
      wshift = 1.0D7 / (conth(ncon) - 109677.576D0/81.0D0**2)
      wmerge = 1.0D7 / (conth(ncon) - emergeh_loc(j))
      IF (wmerge < 0.0D0) wmerge = wshift + wshift
      wcon   = MAX(wshift, wmerge)
      wtail  = 1.0D7 / (1.0D7/wcon - 500.0D0)
      wcon   = MIN(wshift + wshift, wcon)
      IF (wtail < 0.0D0) wtail = wcon + wcon
      wtail  = MIN(wcon + wcon, wtail)
      IF (ifvac == 0) THEN
        wcon  = vac_to_air(wcon)
        wtail = vac_to_air(wtail)
      END IF
      wcon = wcon * dopratio
      IF (ifvac == 0) wl_loc = vac_to_air(1.0D7/(ehyd(nbup)-ehyd(nblo))) * dopratio
      ! Red wing
      IF (wl_loc <= wlend) THEN
        IF (wbegin <= alphahyd(nblo)) THEN
          redcut   = 1.0D7 / (109678.764D0 - 109677.576D0/(DBLE(nbup)-0.8D0)**2 - ehyd(nblo))
          IF (ifvac == 0) redcut = vac_to_air(redcut)
          redcut   = redcut * dopratio
          wlminus1 = 1.0D7 / (ehyd(nbup-1) - ehyd(nblo))
          IF (ifvac == 0) wlminus1 = vac_to_air(wlminus1)
          wlminus1 = wlminus1 * dopratio
          wlminus2 = 1.0D7 / (ehyd(nbup-2) - ehyd(nblo))
          IF (ifvac == 0) wlminus2 = vac_to_air(wlminus2)
          wlminus2 = wlminus2 * dopratio
          kappa0red = kappa0 * REAL( hydrogen_oscillator_strength(nblo,nbup-2)/hydrogen_oscillator_strength(nblo,nbup) / &
                     (ehyd(nbup-2)-ehyd(nblo)) * (ehyd(nbup)-ehyd(nblo)) )
          minred = MAX(1, nbuff)
          wave   = wbegin * ratio**(minred-1)
          DO ibuff = minred, length
            IF (wave >= wcon) THEN
              IF (wave > wlminus1) EXIT
              kappa = kappa0 * hydrogen_line_profile(nblo, nbup, j, wave-wl_loc, dopph)
              IF (wave < wtail) kappa = kappa * REAL((wave-wcon)/(wtail-wcon))
              IF (wave > redcut) THEN
                kappared = kappa0red * hydrogen_line_profile(nblo, nbup-2, j, wave-wlminus2, dopph)
                IF (wave < wtail) kappared = kappared * REAL((wave-wcon)/(wtail-wcon))
                IF (kappared >= kappa) EXIT
              END IF
              buffer(ibuff) = buffer(ibuff) + kappa
              IF (kappa < continuum(ibuff)*cutoff) EXIT
            END IF
            wave = wave * ratio
          END DO
        END IF
        IF (minred == 1) CYCLE line_loop   ! (minred set above)
        IF (wl_loc < wbegin) CYCLE line_loop
      END IF
      ! Blue wing  (613 block)
613   bluecut  = 1.0D7 / (109678.764D0 - 109677.576D0/(DBLE(nbup)+0.8D0)**2 - ehyd(nblo))
      IF (ifvac == 0) bluecut = vac_to_air(bluecut)
      bluecut  = bluecut * dopratio
      wlplus1  = 1.0D7 / (ehyd(nbup+1) - ehyd(nblo))
      IF (ifvac == 0) wlplus1 = vac_to_air(wlplus1)
      wlplus1  = wlplus1 * dopratio
      wlplus2  = 1.0D7 / (ehyd(nbup+2) - ehyd(nblo))
      IF (ifvac == 0) wlplus2 = vac_to_air(wlplus2)
      wlplus2  = wlplus2 * dopratio
      kappa0blue = kappa0 * REAL( hydrogen_oscillator_strength(nblo,nbup+2)/hydrogen_oscillator_strength(nblo,nbup) / &
                   (ehyd(nbup+2)-ehyd(nblo)) * (ehyd(nbup)-ehyd(nblo)) )
      ibuff   = MIN(length+1, nbuff)
      maxblue = ibuff - 1
      wave    = wbegin * ratio**(ibuff-1)
      DO i = 1, maxblue
        ibuff = ibuff - 1
        wave  = wave / ratio
        IF (wave < wcon)    CYCLE line_loop
        IF (wave < wlplus1) CYCLE line_loop
        kappa = kappa0 * hydrogen_line_profile(nblo, nbup, j, wave-wl_loc, dopph)
        IF (wave < wtail) kappa = kappa * REAL((wave-wcon)/(wtail-wcon))
        IF (wave < bluecut) THEN
          kappablue = kappa0blue * hydrogen_line_profile(nblo, nbup+2, j, wave-wlplus2, dopph)
          IF (wave < wtail) kappablue = kappablue * REAL((wave-wcon)/(wtail-wcon))
          IF (kappablue > kappa) CYCLE line_loop
        END IF
        buffer(ibuff) = buffer(ibuff) + kappa
        IF (kappa < continuum(ibuff)*cutoff) CYCLE line_loop
      END DO
      CYCLE line_loop

      ! Alpha/beta isolated red + blue wing  (620/623)
620   IF (wl_loc <= wlend) THEN
        minred = MAX(1, nbuff)
        wave   = wbegin * ratio**(minred-1)
        DO ibuff = minred, length
          kappa = kappa0 * hydrogen_line_profile(nblo, nbup, j, wave-wl_loc, dopph)
          buffer(ibuff) = buffer(ibuff) + kappa
          IF (kappa < continuum(ibuff)*cutoff) EXIT
          wave = wave * ratio
        END DO
        IF (minred == 1) CYCLE line_loop
        IF (wl_loc < wbegin) CYCLE line_loop
      END IF
623   ibuff   = MIN(length+1, nbuff)
      maxblue = ibuff - 1
      wave    = wbegin * ratio**(ibuff-1)
      DO i = 1, maxblue
        ibuff = ibuff - 1
        wave  = wave / ratio
        kappa = kappa0 * hydrogen_line_profile(nblo, nbup, j, wave-wl_loc, dopph)
        buffer(ibuff) = buffer(ibuff) + kappa
        IF (kappa < continuum(ibuff)*cutoff) CYCLE line_loop
      END DO
      CYCLE line_loop

      ! Beta line red wing  (630)
630   IF (wl_loc <= wlend) THEN
        minred = MAX(1, nbuff)
        wave   = wbegin * ratio**(minred-1)
        DO ibuff = minred, length
          IF (wave > alphahyd(nblo)) GOTO 623  ! alpha limit reached -> blue wing
          kappa = kappa0 * hydrogen_line_profile(nblo, nbup, j, wave-wl_loc, dopph)
          buffer(ibuff) = buffer(ibuff) + kappa
          IF (kappa < continuum(ibuff)*cutoff) GOTO 623
          wave = wave * ratio
        END DO
        IF (minred == 1) CYCLE line_loop
        IF (wl_loc < wbegin) CYCLE line_loop
      END IF
      IF (nbuff < 1) CYCLE line_loop
      GOTO 623     ! blue wing

      ! ------------------------------------------------------------------
700   CONTINUE   ! AUTOIONIZING (SHORE) LINE
      ! ------------------------------------------------------------------
      kappa0 = bshore_loc * gf_loc * REAL(xnfpel(nelion))
      kapmin = continuum(MIN(MAX(nbuff,1),length)) * cutoff
      IF (kappa0 < kapmin) CYCLE line_loop
      kappa0 = kappa0 * REAL(EXP(-elo_loc * hckt(j)))
      IF (kappa0 < kapmin) CYCLE line_loop
      IF (linout >= 0) WRITE(15) iline, kappa0
      mlines = mlines + 1
      frelin = REAL(CLIGHT_NM_HZ / wl_loc)
      ! Red wing
      IF (wl_loc <= wlend) THEN
        minred = MAX(1, nbuff)
        freq   = REAL(CLIGHT_NM_HZ / (wbegin * ratio**(minred-1)))
        DO ibuff = minred, length
          epsil = 2.0 * (freq - frelin) / gammar_loc
          kappa = kappa0 * (ashore_loc*epsil + bshore_loc) / (epsil**2 + 1.0) / bshore_loc
          buffer(ibuff) = buffer(ibuff) + kappa
          IF (kappa < continuum(ibuff)*cutoff) EXIT
          freq = freq / REAL(ratio)
        END DO
        IF (nbuff == 1) CYCLE line_loop
        IF (wl_loc < wbegin) CYCLE line_loop
      END IF
      ! Blue wing
      ibuff   = MIN(length+1, nbuff)
      maxblue = ibuff - 1
      freq    = REAL(CLIGHT_NM_HZ / (wbegin * ratio**(ibuff-1)))
      DO i = 1, maxblue
        ibuff = ibuff - 1
        freq  = freq * REAL(ratio)
        epsil = 2.0 * (freq - frelin) / gammar_loc
        kappa = kappa0 * (ashore_loc*epsil + bshore_loc) / (epsil**2 + 1.0) / bshore_loc
        buffer(ibuff) = buffer(ibuff) + kappa
        IF (kappa < continuum(ibuff)*cutoff) CYCLE line_loop
      END DO
      CYCLE line_loop

    END DO line_loop

  END SUBROUTINE compute_line_opacity


  ! ============================================================================
  !  SUBROUTINE run_xnfpelsyn
  !
  !  F90 translation of the XNFPELSYN program (Kurucz).
  !  Replaces the standalone xnfpelsyn executable and the fort.10 file it
  !  produced.  Must be called once before the SYNTHE depth loop.
  !
  !  On entry: mod_atlas_data state is populated by READIN (called by
  !            the host program before calling this subroutine).
  !
  !  On exit:  All xnfpelsyn_shared module arrays are filled:
  !              frqedg_m, wledge_m, cmedge_m, nedge_m
  !              idmol_m, momass_m
  !              freqset_m, numnu_m
  !              xf_t ... xf_xnfh2          (depth structure)
  !              continall_m, contabs_m, contscat_m  (continuum opacities)
  !              xnfpel_m, dopple_m          (ion populations + Doppler widths)
  !              bhyd_gs, bc1_gs, bc2_gs, bsi1_gs, bsi2_gs  (bfudge pops)
  !
  !  ATLAS library routines called (via USE mod_atlas_data):
  !    COMPUTE_ONE_POP -- partition functions and number densities
  !    KAPP            -- continuum opacity (fills ACONT, SIGMAC)
  !    FREEFF          -- free-format number reader
  !
  !  Input:
  !    unit 17  -- edge frequency list (fort.17); read via the ATLAS
  !                free-format reader FREEFF.  Each non-zero value is an
  !                edge frequency, wavelength, or wavenumber (auto-detected
  !                by magnitude); a zero or EOF terminates the list.
  ! ============================================================================
  SUBROUTINE run_xnfpelsyn()

    ! ------------------------------------------------------------------
    !  Access ATLAS state via mod_atlas_data module.
    !  Colliding names imported with _a suffix.
    ! ------------------------------------------------------------------
    USE mod_atlas_data, &
      abund_a => ABUND, ah2p_a => AH2P, ahe1_a => AHE1, ahe2_a => AHE2, &
      ahemin_a => AHEMIN, ahline_a => AHLINE, ahmin_a => AHMIN, &
      ahyd_a => AHYD, atmass_a => ATMASS, bhyd_a => BHYD, bmin_a => BMIN, &
      elem_a => ELEM, hckt_a => HCKT, hkt_a => HKT, ifturb_a => IFTURB, &
      itemp_a => ITEMP, nrhox_a => NRHOX, p_a => P, ptotal_a => PTOTAL, &
      pturb_a => PTURB, rho_a => RHO, rhox_a => RHOX, shline_a => SHLINE, &
      shmin_a => SHMIN, shyd_a => SHYD, sigh2_a => SIGH2, sigh_a => SIGH, &
      sighe_a => SIGHE, t_a => T, tk_a => TK, tkev_a => TKEV, tlog_a => TLOG, &
      trbcon_a => TRBCON, trbfdg_a => TRBFDG, trbpow_a => TRBPOW, &
      trbsnd_a => TRBSND, vturb_a => VTURB, xnatom_a => XNATOM, xne_a => XNE, &
      XNF_a => XNF, XNFP_a => XNFP

    ! Population output arrays (local; filled by COMPUTE_ONE_POP)
    REAL(8) :: xnfh_x(kw), xnfhe_x(kw,2)
    REAL(8) :: xnfph_x(kw,2), xnfphe_x(kw,3)
    REAL(8) :: xnfpb_x(kw,1), xnfpc_x(kw,2), xnfpo_x(kw,1)
    REAL(8) :: xnfpna_x(kw,1), xnfpmg_x(kw,2), xnfpal_x(kw,2)
    REAL(8) :: xnfpsi_x(kw,2), xnfpk_x(kw,1), xnfpca_x(kw,2)
    REAL(8) :: xnfpfe_x(kw,1)
    REAL(8) :: xnfpch_x(kw), xnfpoh_x(kw)

    !  Local variables
    ! ------------------------------------------------------------------
    INTEGER  :: i, j, nu, in_e, last1, nelem, ion, nsteps_x
    REAL(8)  :: edge, sav, freq15_x, wave_x, waveno_lc, stepwt_x
    REAL(8)  :: eq, rco
    REAL(8)  :: a_sort(me)              ! sortable wavelength absolute values
    REAL(8)  :: xnfh2_lc(kw)           ! local H2 number density
    REAL(8)  :: xnfph2_lc(kw)          ! local H2 proton-fraction density
    REAL(8)  :: xnfpco_lc(kw)          ! local CO number density
    REAL(8)  :: xnfp_lc(kw, 10, mw)    ! full ion population array
    REAL(8)  :: xnfpel_lc(6, mw)       ! per-depth slice for writing
    REAL(8)  :: dopple_lc(6, mw)       ! per-depth Doppler widths
    CHARACTER(LEN=1) :: card_x(81)     ! input line buffer for FREEFF

    ! ------------------------------------------------------------------
    !  IDMOL and MOMASS data (molecule codes and masses)
    !  These reproduce the DATA statements from xnfpelsyn.for exactly.
    ! ------------------------------------------------------------------
    REAL(8), PARAMETER :: idmol_data(mm) = (/ &
       101D0,  106D0,  107D0,  108D0,  606D0,  607D0,  608D0,  707D0,  708D0, &
       808D0,  112D0,  113D0,  114D0,  812D0,  813D0,  814D0,  116D0,  120D0, &
       816D0,  820D0,  821D0,  822D0,  823D0,  103D0,  104D0,  105D0,  109D0, &
       115D0,  117D0,  121D0,  122D0,  123D0,  124D0,  125D0,  126D0, 106.01D0, &
    107.01D0, 108.01D0, 112.01D0, 113.01D0, 114.01D0, 120.01D0, 111D0,  119D0, 10101.01D0, &
       817D0,  824D0,  825D0,  826D0, 10108D0, 60808D0, 10106D0, 60606D0,  127D0, &
       128D0,  129D0,  827D0,  828D0,  829D0, 608.01D0,  408D0,  508D0,  815D0, &
     10808D0, 10811D0, 10812D0, 10820D0, 10106D0, 10107D0, 10116D0, 10606D0, 10607D0, &
     10608D0, 10708D0, 60816D0, 61616D0, 70708D0, 70808D0, 80814D0, 80816D0, 1010106D0, &
    1010107D0, 1010606D0, 101010106D0, 101010114D0, 614D0, 60614D0, 60607D0, &
    6060707D0, 6060607D0, &
       839D0,  840D0,  857D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0 /)

    REAL(8), PARAMETER :: momass_data(mm) = (/ &
       2D0,  13D0,  15D0,  17D0,  24D0,  26D0,  28D0,  28D0,  30D0, &
      32D0,  25D0,  28D0,  29D0,  40D0,  43D0,  44D0,  33D0,  41D0, &
      48D0,  56D0,  61D0,  64D0,  67D0,   8D0,  10D0,  12D0,  20D0, &
      32D0,  36D0,  46D0,  49D0,  52D0,  53D0,  56D0,  57D0,  13D0, &
      15D0,  17D0,  25D0,  28D0,  29D0,  41D0,  24D0,  40D0,   3D0, &
      51D0,  68D0,  71D0,  72D0,  18D0,  44D0,  14D0,  36D0,  60D0, &
      59D0,  64D0,  75D0,  74D0,  79D0,  28D0,  25D0,  27D0,  47D0, &
      33D0,  30D0,  31D0,  57D0,  14D0,  16D0,  34D0,  25D0,  27D0, &
      29D0,  31D0,  60D0,  74D0,  44D0,  46D0,  60D0,  64D0,  15D0, &
      17D0,  26D0,  16D0,  32D0,  40D0,  52D0,  38D0,  52D0,  50D0, &
     104D0, 107D0, 155D0, 999D0, 999D0, 999D0, 999D0, 999D0, 999D0, 999D0 /)

    ! ------------------------------------------------------------------
    !  Initialise
    ! ------------------------------------------------------------------
    idmol_m  = idmol_data
    momass_m = momass_data
    card_x   = ' '

    ! ------------------------------------------------------------------
    !  ATLAS setup flags  (mirror of the original program's flag setting)
    ! ------------------------------------------------------------------
    IFOP(14) = 0
    IFOP(15) = 0
    IFOP(16) = 0
    IFOP(17) = 0
    MORE     = 1
    MAXPOW   = 99
    LAST     = 81
    NUMCOL   = 1
    itemp_a    = 0
    ! Respect the IFMOL flag set by READIN from the model control cards.
    ! When IFMOL=1 ("MOLECULES ON"), we need IFPRES=1 so that NMOLEC
    ! runs inside COMPUTE_ONE_POP and populates the molecular number
    ! density arrays.  When IFMOL=0, IFPRES=0 is fine: COMPUTE_ONE_POP
    ! takes the atomic-only path via PFSAHA.
    IF (IFMOL == 1) THEN
      IFPRES = 1
    ELSE
      IFPRES = 0
    END IF
    ITER     = 1
    NUMITS   = 1

    ! ------------------------------------------------------------------
    !  SECTION 1.  READ EDGE FREQUENCY GRID FROM UNIT 17
    !
    !  Each call to XFREEF (ATLAS free-format reader) returns one numeric
    !  value from fort.17.  Values are auto-detected as wavelengths (nm),
    !  wavenumbers (cm^{-1}), or frequencies (Hz) by magnitude:
    !    |value| < 1e6    -> wavelength in nm
    !    1e6 <= |v| < 1e25 -> frequency in Hz
    !    |value| >= 1e25  -> wavenumber in cm^{-1} (multiply by 1e25 first)
    !  A value of exactly 0 ends the list.
    ! ------------------------------------------------------------------
    in_e = 0
    edge_read: DO
      edge = xfreef_x(card_x)
      IF (edge == 0.0D0) EXIT edge_read
      in_e = in_e + 1
      IF (in_e > me) THEN
        WRITE(6,'(A,I5)') ' run_xnfpelsyn: too many edges, max=', me
        STOP
      END IF

      IF (ABS(edge) >= 1.0D25) THEN
        ! Wavenumber (cm^{-1}): convert to nm and Hz
        cmedge_m(in_e) = edge / 1.0D25
        wledge_m(in_e) = 1.0D7 / cmedge_m(in_e)
        frqedg_m(in_e) = CLIGHT_NM_HZ / wledge_m(in_e)
      ELSE IF (ABS(edge) < 1.0D6) THEN
        ! Wavelength in nm
        wledge_m(in_e) = edge
        cmedge_m(in_e) = 1.0D7 / edge
        frqedg_m(in_e) = CLIGHT_NM_HZ / wledge_m(in_e)
      ELSE
        ! Frequency in Hz
        frqedg_m(in_e) = edge
        wledge_m(in_e) = CLIGHT_NM_HZ / edge
        cmedge_m(in_e) = 1.0D7 / wledge_m(in_e)
      END IF
      a_sort(in_e) = ABS(wledge_m(in_e))
    END DO edge_read

    ! Sort edges by ascending wavelength (bubble sort; list is short)
    DO i = 2, in_e
      last1 = in_e - i + 2
      DO j = 2, last1
        IF (a_sort(j) >= a_sort(j-1)) CYCLE
        sav = a_sort(j-1);    a_sort(j-1)    = a_sort(j);    a_sort(j)    = sav
        sav = frqedg_m(j-1);  frqedg_m(j-1)  = frqedg_m(j);  frqedg_m(j)  = sav
        sav = wledge_m(j-1);  wledge_m(j-1)  = wledge_m(j);  wledge_m(j)  = sav
        sav = cmedge_m(j-1);  cmedge_m(j-1)  = cmedge_m(j);  cmedge_m(j)  = sav
      END DO
    END DO
    nedge_m = in_e

    ! Build the 3-point continuum frequency grid between edges
    numnu_m = 0
    DO i = 1, nedge_m - 1
      numnu_m = numnu_m + 1
      freqset_m(numnu_m) = ABS(frqedg_m(i))   / 1.0000001D0
      numnu_m = numnu_m + 1
      freqset_m(numnu_m) = CLIGHT_NM_HZ / (ABS(wledge_m(i)) + ABS(wledge_m(i+1))) * 2.0D0
      numnu_m = numnu_m + 1
      freqset_m(numnu_m) = ABS(frqedg_m(i+1)) * 1.0000001D0
    END DO

    ! ------------------------------------------------------------------
    !  SECTION 2.  POPULATE ALL ION/MOLECULAR NUMBER DENSITIES
    !
    !  COMPUTE_ALL_POPS fills the mod_atlas_data arrays XNF and XNFP
    !  for all elements Z=1-99, H2, CO, and (if IFMOL=1) molecules.
    !  KAPP's opacity routines (HOP, HMINOP, etc.) read directly from
    !  XNF/XNFP, so this single call provides everything KAPP needs.
    ! ------------------------------------------------------------------
    itemp_a = itemp_a + 1

    ! When IFMOL=1, NMOLEC must be called with MODE=1 (or 11) to
    ! populate the XNFPMOL array (number density / partition function).
    ! COMPUTE_ALL_POPS triggers NMOLEC via COMPUTE_ONE_POP, but its
    ! first call uses MODE=12, which tells NMOLEC to skip the XNFPMOL
    ! computation.  Subsequent calls find ITEMP_PREV==ITEMP and don't
    ! re-trigger NMOLEC at all.  So we call NMOLEC(1) explicitly here.
    !
    ! NMOLEC also recalculates XNE as part of the molecular + ionisation
    ! equilibrium.  The result is self-consistent but can differ slightly
    ! (<1%) from the converged XNE stored in the model atmosphere.  We
    ! save and restore XNE to keep the model's converged values.
    IF (IFMOL == 1) THEN
      BLOCK
        REAL(8) :: xne_save(kw)
        xne_save(1:nrhox_a) = xne_a(1:nrhox_a)
        CALL NMOLEC(1)
        xne_a(1:nrhox_a) = xne_save(1:nrhox_a)
      END BLOCK
    END IF

    CALL COMPUTE_ALL_POPS()

    ! Extract specific populations needed by later sections.
    ! NOFF(IZ) = IZ*(IZ+1)/2: H=1, He=3, C=21, Si=105
    ! XNFP_a(j, NOFF(IZ)+ion-1) = mode-11 result for element IZ, stage ion
    ! XNF_a(j, NOFF(IZ)+ion-1)  = mode-12 result
    DO j = 1, nrhox_a
      xnfh_x(j)      = XNF_a(j, 1)          ! H mode-12
      xnfhe_x(j, 1)  = XNF_a(j, 3)          ! He I mode-12
      xnfhe_x(j, 2)  = XNF_a(j, 4)          ! He II mode-12
      xnfph_x(j, 1)  = XNFP_a(j, 1)         ! H I mode-11
      xnfph_x(j, 2)  = XNFP_a(j, 2)         ! H II mode-11
      xnfphe_x(j, 1) = XNFP_a(j, 3)         ! He I mode-11
      xnfphe_x(j, 2) = XNFP_a(j, 4)         ! He II mode-11
      xnfphe_x(j, 3) = XNFP_a(j, 5)         ! He III mode-11
      xnfpc_x(j, 1)  = XNFP_a(j, 21)        ! C I mode-11
      xnfpc_x(j, 2)  = XNFP_a(j, 22)        ! C II mode-11
      xnfpsi_x(j, 1) = XNFP_a(j, 105)       ! Si I mode-11
      xnfpsi_x(j, 2) = XNFP_a(j, 106)       ! Si II mode-11
      xnfpo_x(j, 1)  = XNFP_a(j, 36)        ! O I mode-11
    END DO

    ! ------------------------------------------------------------------
    !  SECTION 2b.  GROUND-STATE BOLTZMANN POPULATIONS FOR BFUDGE
    !
    !  The source-function fudge factor bfudge_sv requires ground-state
    !  Boltzmann populations: bhyd(j,1), bc1(j,1)/bc2(j,1), bsi1(j,1)/bsi2(j,1).
    !  In the old F77 code these were stored in COMMON B-arrays as a side effect
    !  of PFSAHA.  The modernized PFSAHA returns only ion-level totals via its
    !  argument, so we compute the ground-state values explicitly here using
    !  the PFSAHA level energies and degeneracies.
    !
    !  For each species, the ground-state population is:
    !    B_gs(j) = g_1 * exp(-E_1 * HCKT(j)) / U(ion) * F(ion) * XNATOM(j) * XABUND(j,IZ)
    !  where the mode-11 COMPUTE_ONE_POP result already gives:
    !    xnfp(j,ion) = F(ion) / U(ion) * XNATOM(j) * XABUND(j,IZ)
    !  so:
    !    B_gs(j) = g_1 * exp(-E_1 * HCKT(j)) * xnfp(j,ion)
    !
    !  Level data from PFSAHA (cm^{-1} and statistical weights):
    !    H I:   E=0,       g=2     (EHYD(1), GHYD(1))
    !    C I:   E=29.60,   g=9     (EC1(1),  GC1(1))
    !    C II:  E=42.48,   g=6     (EC2(1),  GC2(1))
    !    Si I:  E=149.681, g=9     (ESI1(1), GSI1(1))
    !    Si II: E=191.55,  g=6     (ESI2(1), GSI2(1))
    ! ------------------------------------------------------------------
    DO j = 1, nrhox_a
      bhyd_gs(j)  = 2.0D0 * EXP(-0.0D0    * hckt_a(j)) * xnfph_x(j,1)
      bc1_gs(j)   = 9.0D0 * EXP(-29.60D0  * hckt_a(j)) * xnfpc_x(j,1)
      bc2_gs(j)   = 6.0D0 * EXP(-42.48D0  * hckt_a(j)) * xnfpc_x(j,2)
      bsi1_gs(j)  = 9.0D0 * EXP(-149.681D0* hckt_a(j)) * xnfpsi_x(j,1)
      bsi2_gs(j)  = 6.0D0 * EXP(-191.55D0 * hckt_a(j)) * xnfpsi_x(j,2)
    END DO

    ! ------------------------------------------------------------------
    !  SECTION 3.  H2 EQUILIBRIUM POPULATIONS
    !
    !  Polynomial fit to the H2 partition function (valid T < 9000 K).
    !  Coefficients reproduce the Kurucz fit exactly.
    ! ------------------------------------------------------------------
    xnfh2_lc  = 0.0D0
    xnfph2_lc = 0.0D0
    xnfpco_lc = 0.0D0
    DO j = 1, nrhox_a
      IF (t_a(j) > 9000.0D0) CYCLE
      eq = EXP(4.478D0/tkev_a(j) - 4.64584D1 + &
          (1.63660D-3 + (-4.93992D-7 + (1.11822D-10 + (-1.49567D-14 + &
          (1.06206D-18 - 3.08720D-23*t_a(j))*t_a(j))*t_a(j))*t_a(j))*t_a(j))*t_a(j) - &
           1.5D0*tlog_a(j))
      xnfh2_lc(j)  = xnfh_x(j)**2 * eq
      xnfph2_lc(j) = xnfph_x(j,1)**2 * eq
      xnfpco_lc(j) = xnfpc_x(j,1) * xnfpo_x(j,1) * &
          EXP(11.091D0/tkev_a(j) - 49.0414D0 + &
           14.0306D-4*t_a(j) - 26.6341D-8*t_a(j)**2 + &
           35.382D-12*t_a(j)**3 - 26.5424D-16*t_a(j)**4 + &
           8.32385D-20*t_a(j)**5 - 1.5D0*tlog_a(j))
    END DO

    ! ------------------------------------------------------------------
    !  SECTION 4.  CONTINUUM OPACITY LOOP
    !
    !  For each of the NUMNU_M grid points, compute the continuum opacity
    !  at all depths using KAPP (ATLAS library), then store the result.
    ! ------------------------------------------------------------------
    DO nu = 1, numnu_m
      FREQ   = freqset_m(nu)
      freq15_x = FREQ / 1.0D15
      FREQLG = LOG(FREQ)
      rco      = 0.0D0
      WAVE   = CLIGHT_NM_HZ / FREQ
      WAVENO = 1.0D7 / WAVE
      DO j = 1, nrhox_a
        EHVKT(j) = EXP(-FREQ * hkt_a(j))
        STIM(j)  = 1.0D0 - EHVKT(j)
        BNU(j)   = 1.47439D-2 * freq15_x**3 * EHVKT(j) / STIM(j)
      END DO
      CALL kapp()
      DO j = 1, nrhox_a
        continall_m(nu,j) = LOG10(ACONT(j) + SIGMAC(j))
        contabs_m  (nu,j) = LOG10(ACONT(j))
        contscat_m (nu,j) = LOG10(SIGMAC(j))
      END DO
    END DO

    ! ------------------------------------------------------------------
    !  SECTION 5.  STORE DEPTH STRUCTURE
    !
    !  Copy the ATLAS COMMON arrays into the module depth-structure arrays.
    !  These replace the single depth-structure record in fort.10.
    ! ------------------------------------------------------------------
    DO j = 1, nrhox_a
      xf_t(j)       = t_a(j)
      xf_tkev(j)    = tkev_a(j)
      xf_tk(j)      = tk_a(j)
      xf_hkt(j)     = hkt_a(j)
      xf_tlog(j)    = tlog_a(j)
      xf_hckt(j)    = hckt_a(j)
      xf_p(j)       = p_a(j)
      xf_xne(j)     = xne_a(j)
      xf_xnatom(j)  = xnatom_a(j)
      xf_rho(j)     = rho_a(j)
      xf_rhox(j)    = rhox_a(j)
      xf_vturb(j)   = vturb_a(j)
      xf_xnfh(j)    = xnfh_x(j)
      xf_xnfhe(j,1) = xnfhe_x(j,1)
      xf_xnfhe(j,2) = xnfhe_x(j,2)
      xf_xnfh2(j)   = xnfh2_lc(j)
    END DO

    ! ------------------------------------------------------------------
    !  SECTION 6.  FULL ION POPULATION ARRAY  (XNFP)
    !
    !  Build the complete (kw, 10, mw) array of number densities for all
    !  elements (Z=1..99) and all ionisation stages, plus molecules.
    !  This mirrors the long sequence of POPS calls in the original.
    ! ------------------------------------------------------------------
    xnfp_lc = 0.0D0

    ! Elements Z=1..28 with explicit ion-stage codes
    CALL COMPUTE_ONE_POP(1.01D0,  11, xnfp_lc(1,1,1))
    CALL COMPUTE_ONE_POP(2.02D0,  11, xnfp_lc(1,1,2))
    CALL COMPUTE_ONE_POP(3.03D0,  11, xnfp_lc(1,1,3))
    CALL COMPUTE_ONE_POP(4.03D0,  11, xnfp_lc(1,1,4))
    CALL COMPUTE_ONE_POP(5.03D0,  11, xnfp_lc(1,1,5))
    CALL COMPUTE_ONE_POP(6.05D0,  11, xnfp_lc(1,1,6))
    CALL COMPUTE_ONE_POP(7.05D0,  11, xnfp_lc(1,1,7))
    CALL COMPUTE_ONE_POP(8.05D0,  11, xnfp_lc(1,1,8))
    CALL COMPUTE_ONE_POP(9.05D0,  11, xnfp_lc(1,1,9))
    CALL COMPUTE_ONE_POP(10.05D0, 11, xnfp_lc(1,1,10))
    CALL COMPUTE_ONE_POP(11.05D0, 11, xnfp_lc(1,1,11))
    CALL COMPUTE_ONE_POP(12.05D0, 11, xnfp_lc(1,1,12))
    CALL COMPUTE_ONE_POP(13.05D0, 11, xnfp_lc(1,1,13))
    CALL COMPUTE_ONE_POP(14.05D0, 11, xnfp_lc(1,1,14))
    CALL COMPUTE_ONE_POP(15.05D0, 11, xnfp_lc(1,1,15))
    CALL COMPUTE_ONE_POP(16.05D0, 11, xnfp_lc(1,1,16))
    CALL COMPUTE_ONE_POP(17.04D0, 11, xnfp_lc(1,1,17))
    CALL COMPUTE_ONE_POP(18.04D0, 11, xnfp_lc(1,1,18))
    CALL COMPUTE_ONE_POP(19.04D0, 11, xnfp_lc(1,1,19))
    CALL COMPUTE_ONE_POP(20.09D0, 11, xnfp_lc(1,1,20))
    CALL COMPUTE_ONE_POP(21.09D0, 11, xnfp_lc(1,1,21))
    CALL COMPUTE_ONE_POP(22.09D0, 11, xnfp_lc(1,1,22))
    CALL COMPUTE_ONE_POP(23.09D0, 11, xnfp_lc(1,1,23))
    CALL COMPUTE_ONE_POP(24.09D0, 11, xnfp_lc(1,1,24))
    CALL COMPUTE_ONE_POP(25.09D0, 11, xnfp_lc(1,1,25))
    CALL COMPUTE_ONE_POP(26.09D0, 11, xnfp_lc(1,1,26))
    CALL COMPUTE_ONE_POP(27.09D0, 11, xnfp_lc(1,1,27))
    CALL COMPUTE_ONE_POP(28.09D0, 11, xnfp_lc(1,1,28))

    ! Elements Z=29..99
    DO nelem = 29, 99
      CALL COMPUTE_ONE_POP(DBLE(nelem) + 0.02D0, 11, xnfp_lc(1,1,nelem))
    END DO

    ! H2 and CO special entries
    DO j = 1, nrhox_a
      xnfp_lc(j, 6, 40) = xnfph2_lc(j)
      xnfp_lc(j, 6, 46) = xnfpco_lc(j)
    END DO

    ! Molecules (NELEM=40..mw)
    IF (IFMOL == 1) THEN
      DO nelem = 40, mw
        CALL COMPUTE_ONE_POP(idmol_data(nelem-39), 1, xnfp_lc(1,6,nelem))
      END DO
    END IF

    ! Higher ionisation stages for elements 20..28 (stages 7-10 -> 5,30+N etc.)
    DO j = 1, nrhox_a
      DO nelem = 20, 28
        xnfp_lc(j, 5, 30+nelem) = xnfp_lc(j, 7, nelem)
        xnfp_lc(j, 5, 40+nelem) = xnfp_lc(j, 8, nelem)
        xnfp_lc(j, 5, 50+nelem) = xnfp_lc(j, 9, nelem)
        xnfp_lc(j, 5, 60+nelem) = xnfp_lc(j, 10, nelem)
      END DO
    END DO

    ! ------------------------------------------------------------------
    !  SECTION 7.  PACK XNFPEL / DOPPLE INTO MODULE ARRAYS
    !
    !  Compute Doppler widths and store both populations and Doppler
    !  widths into the per-depth module arrays xnfpel_m / dopple_m.
    ! ------------------------------------------------------------------
    DO j = 1, nrhox_a

      DO nelem = 1, mw
        DO ion = 1, 6
          xnfpel_lc(ion, nelem) = xnfp_lc(j, ion, nelem)
        END DO
      END DO

      ! Atomic Doppler width: sqrt(2kT/m + v_turb^2) / c
      ! TK = k_B * T (ergs); ATMASS in amu; 1 amu = 1.660e-24 g
      DO nelem = 1, 99
        dopple_lc(1, nelem) = SQRT(2.0D0 * tk_a(j) / &
                              (atmass_a(nelem) * 1.660D-24) + &
                               vturb_a(j)**2) / 2.99792458D10
        DO ion = 2, 6
          dopple_lc(ion, nelem) = dopple_lc(1, nelem)
        END DO
      END DO

      ! High-ionisation stage aliases (same thermal width as parent element)
      DO nelem = 20, 28
        dopple_lc(5, 30+nelem) = dopple_lc(1, nelem)
        dopple_lc(5, 40+nelem) = dopple_lc(1, nelem)
        dopple_lc(5, 50+nelem) = dopple_lc(1, nelem)
        dopple_lc(5, 60+nelem) = dopple_lc(1, nelem)
      END DO

      ! Molecular Doppler widths (use MOMASS for reduced mass)
      DO nelem = 40, mw
        dopple_lc(6, nelem) = SQRT(2.0D0 * tk_a(j) / &
                              (momass_data(nelem-39) * 1.660D-24) + &
                               vturb_a(j)**2) / 2.99792458D10
      END DO

      ! Store into module arrays
      DO nelem = 1, mw
        DO ion = 1, 6
          xnfpel_m(ion, nelem, j) = xnfpel_lc(ion, nelem)
          dopple_m (ion, nelem, j) = dopple_lc(ion, nelem)
        END DO
      END DO

    END DO   ! depth loop

  CONTAINS

    ! ---------------------------------------------------------------
    !  XFREEF_X  -- thin wrapper around the ATLAS FREEFF free-format
    !  reader.  When FREEFF signals that a new input line is needed
    !  (IFFAIL /= 0), read the next line from unit 17 and retry.
    ! ---------------------------------------------------------------
    FUNCTION xfreef_x(card) RESULT(val)
      CHARACTER(LEN=1), INTENT(INOUT) :: card(81)
      REAL(8) :: val
      INTEGER :: l

      MORE  = 1
      val     = freeff(card)
      IF (IFFAIL == 0) RETURN
      l = LAST - 1
      READ(17, '(80A1)', END=99) (card(i), i=1,l)
      NUMCOL = 1
      val      = freeff(card)
      RETURN
99    val = 0.0D0   ! EOF -> signal end of edge list
    END FUNCTION xfreef_x

  END SUBROUTINE run_xnfpelsyn


  ! ==========================================================================
END MODULE synthe_module
