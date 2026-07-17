! joshtest: standalone texture probe for the JOSH-path numerics
! (MAP1 remap, DERIV slopes, INTEG/PARCOE quadrature), following the
! pftest pattern of calling the ACTUAL production routines.
!
! Reads josh_profile.txt (N; then RHOX, T, ABROSS per layer) extracted
! from a converged model.  Each test sweeps a perturbation amplitude
! eps over an analytic (infinitely differentiable) family, so any
! non-smoothness of the recorded outputs vs eps is method noise, not
! signal.  Results go to josh_fortran.csv for josh_verify.py.
!
!   A  : T -> T*(1+eps*g(tau)), S=Planck(T); MAP1 of S onto the fixed
!        51-point XTAU8 grid; records surface flux H = sum(CH_J*XS8)
!        and the raw interpolant at tau=0.5.  Grids fixed, only data
!        moves: isolates the data-dependence of MAP1's blend weights.
!   A2 : opacity -> ABROSS*(1+eps*g), S fixed; INTEG rebuilds TAUNU so
!        the source grid slides under XTAU8: interval-crossing noise.
!   B  : same S(eps) as A; HNU = DERIV(S)/3, JMINS = DERIV(HNU) --
!        the deep-atmosphere Eddington branch's nested derivative.
!   C  : accuracy of INTEG against an analytic integral exp(a*x) on
!        the real RHOX grid (no sweep).
PROGRAM JOSHTEST
  USE mod_atlas_data, ONLY: MAP1, DERIV, INTEG, IQUAD
  IMPLICIT NONE

  INTEGER, PARAMETER :: NXTAU = 51, NEPS = 401
  REAL(8), PARAMETER :: EPSMAX = 2.0D-3
  REAL(8), PARAMETER :: HPL = 6.62607015D-27, KB = 1.380649D-16, &
                        CL = 2.99792458D10
  REAL(8), PARAMETER :: WAVECM = 500.0D-7   ! 500 nm probe frequency

  ! XTAU8 and CH_J duplicated verbatim from JOSH (local PARAMETERs there)
  REAL(8), PARAMETER :: XTAU8(NXTAU) = (/ &
    0.0D0, 0.0000032D0, 0.0000056D0, 0.00001D0, 0.000018D0, &
    0.000032D0, 0.000056D0, 0.0001D0, 0.00018D0, 0.00032D0, &
    0.00056D0, 0.001D0, 0.0018D0, 0.0032D0, 0.0056D0, &
    0.01D0, 0.016D0, 0.025D0, 0.042D0, 0.065D0, &
    0.096D0, 0.139D0, 0.196D0, 0.273D0, 0.375D0, &
    0.5D0, 0.63D0, 0.78D0, 0.95D0, 1.15D0, &
    1.35D0, 1.6D0, 1.85D0, 2.15D0, 2.45D0, &
    2.75D0, 3.15D0, 3.65D0, 4.25D0, 5.0D0, &
    6.0D0, 7.0D0, 8.0D0, 9.0D0, 10.0D0, &
    11.5D0, 13.0D0, 14.5D0, 16.0D0, 18.0D0, 20.0D0 /)
  REAL(8), PARAMETER :: CH_J(NXTAU) = (/ &
    7.15528131D-07, 1.49142693D-06, 1.52106577D-06, 2.98150826D-06, &
    5.33941056D-06, 9.13329677D-06, 1.61715943D-05, 2.97035986D-05, &
    5.33166603D-05, 9.11154202D-05, 1.61084638D-04, 2.95118050D-04, &
    5.27450291D-04, 8.67939554D-04, 1.61498412D-03, 2.50720908D-03, &
    3.20994272D-03, 5.61912498D-03, 8.60872678D-03, 1.04706492D-02, &
    1.33110350D-02, 1.62635669D-02, 1.90288834D-02, 2.18877215D-02, &
    2.36015432D-02, 2.10819542D-02, 1.80345085D-02, 1.64786074D-02, &
    1.49382707D-02, 1.19676525D-02, 9.90213640D-03, 8.17766134D-03, &
    6.11252524D-03, 4.84035723D-03, 3.06078210D-03, 2.40512565D-03, &
    2.01712688D-03, 1.33288081D-03, 7.83530239D-04, 4.31428343D-04, &
    1.76504589D-04, 4.75738016D-05, 1.65963702D-05, 5.41117970D-06, &
    2.08043571D-06, 7.11612643D-07, 8.08788982D-08, 1.95130507D-08, &
    4.33638281D-09, 8.87765583D-10, 3.90236420D-11 /)

  INTEGER :: N, J, IE, IDUM, LUIN, LUOUT, JP(3), K
  REAL(8) :: EPS, H, FREQ, AEXP, TARGETS(3), DBEST
  REAL(8), ALLOCATABLE :: RHOX(:), T(:), ABR(:), TAUROS(:), G(:), &
       SNU(:), SNU0(:), TP(:), TAUNU(:), ABP(:), HNU(:), JM(:), &
       F(:), FI(:)
  REAL(8) :: XS8(NXTAU)

  FREQ = CL / WAVECM

  OPEN(NEWUNIT=LUIN, FILE='josh_profile.txt', STATUS='OLD', ACTION='READ')
  READ(LUIN, *) N
  ALLOCATE(RHOX(N), T(N), ABR(N), TAUROS(N), G(N), SNU(N), SNU0(N), &
           TP(N), TAUNU(N), ABP(N), HNU(N), JM(N), F(N), FI(N))
  DO J = 1, N
    READ(LUIN, *) RHOX(J), T(J), ABR(J)
  END DO
  CLOSE(LUIN)

  ! Rosseland tau scale exactly as production (ROSS mode 3)
  CALL INTEG(RHOX, ABR, TAUROS, N, ABR(1)*RHOX(1))

  ! Smooth Gaussian bump in log10(tau), centered on tau=1, width 0.35 dex
  DO J = 1, N
    G(J) = exp(-0.5D0 * ((log10(TAUROS(J)) - 0.0D0) / 0.35D0)**2)
  END DO

  ! Probe layers nearest tau = 0.5, 2, 10
  TARGETS = (/ 0.5D0, 2.0D0, 10.0D0 /)
  DO K = 1, 3
    JP(K) = 1
    DBEST = huge(1.0D0)
    DO J = 1, N
      IF (abs(TAUROS(J) - TARGETS(K)) .LT. DBEST) THEN
        DBEST = abs(TAUROS(J) - TARGETS(K))
        JP(K) = J
      END IF
    END DO
  END DO

  OPEN(NEWUNIT=LUOUT, FILE='josh_fortran.csv', STATUS='REPLACE', ACTION='WRITE')
  WRITE(LUOUT, '(A,3(",",I4))') 'PROBES', JP(1), JP(2), JP(3)
  DO J = 1, NXTAU
    WRITE(LUOUT, '(A,",",I4,2(",",ES24.16))') 'G', J, XTAU8(J), CH_J(J)
  END DO
  DO J = 1, N
    WRITE(LUOUT, '(A,",",I4,4(",",ES24.16))') 'P', J, RHOX(J), T(J), ABR(J), TAUROS(J)
  END DO

  !-------------------------------------------------------------------
  ! Test A: fixed grids, S(eps) analytic; MAP1 onto XTAU8 + CH_J flux
  !-------------------------------------------------------------------
  DO IE = 0, NEPS - 1
    EPS = EPSMAX * IE / DBLE(NEPS - 1)
    DO J = 1, N
      TP(J)  = T(J) * (1.0D0 + EPS * G(J))
      SNU(J) = BNUF(TP(J))
    END DO
    XS8 = 0.0D0
    IDUM = MAP1(TAUROS, SNU, N, XTAU8, XS8, NXTAU)
    H = SUM(CH_J * XS8)
    WRITE(LUOUT, '(A,3(",",ES24.16))') 'A', EPS, H, XS8(26)
  END DO

  !-------------------------------------------------------------------
  ! Test A2: S fixed, opacity(eps) analytic; INTEG rebuilds TAUNU
  !-------------------------------------------------------------------
  DO J = 1, N
    SNU0(J) = BNUF(T(J))
  END DO
  DO IE = 0, NEPS - 1
    EPS = EPSMAX * IE / DBLE(NEPS - 1)
    ABP = ABR * (1.0D0 + EPS * G)
    CALL INTEG(RHOX, ABP, TAUNU, N, ABP(1)*RHOX(1))
    XS8 = 0.0D0
    IDUM = MAP1(TAUNU, SNU0, N, XTAU8, XS8, NXTAU)
    H = SUM(CH_J * XS8)
    WRITE(LUOUT, '(A,5(",",ES24.16))') 'A2', EPS, H, &
      TAUNU(JP(1)), TAUNU(JP(2)), TAUNU(JP(3))
  END DO

  !-------------------------------------------------------------------
  ! Test B: nested DERIV (deep-branch Eddington flux and J-S)
  !-------------------------------------------------------------------
  DO IE = 0, NEPS - 1
    EPS = EPSMAX * IE / DBLE(NEPS - 1)
    DO J = 1, N
      TP(J)  = T(J) * (1.0D0 + EPS * G(J))
      SNU(J) = BNUF(TP(J))
    END DO
    CALL DERIV(TAUROS, SNU, HNU, N)
    HNU = HNU / 3.0D0
    CALL DERIV(TAUROS, HNU, JM, N)
    WRITE(LUOUT, '(A,7(",",ES24.16))') 'B', EPS, &
      HNU(JP(1)), HNU(JP(2)), HNU(JP(3)), JM(JP(1)), JM(JP(2)), JM(JP(3))
  END DO

  !-------------------------------------------------------------------
  ! Test C: INTEG accuracy vs analytic integral of exp(a*x).
  ! Exponent normalized so F spans 5 decades on any model's RHOX grid.
  !-------------------------------------------------------------------
  AEXP = log(1.0D5) / (RHOX(N) - RHOX(1))
  WRITE(LUOUT, '(A,",",ES24.16)') 'CA', AEXP
  F = exp(AEXP * RHOX)
  CALL INTEG(RHOX, F, FI, N, 0.0D0)
  DO J = 1, N
    WRITE(LUOUT, '(A,",",I4,3(",",ES24.16))') 'C', J, RHOX(J), &
      (exp(AEXP*RHOX(J)) - exp(AEXP*RHOX(1))) / AEXP, FI(J)
  END DO

  ! Same integral through the production INTEG dispatch with quad=1
  ! (Steffen path), for cross-checking against the independent Python
  ! implementation in josh_verify.py.
  IQUAD = 1
  CALL INTEG(RHOX, F, FI, N, 0.0D0)
  IQUAD = 0
  DO J = 1, N
    WRITE(LUOUT, '(A,",",I4,3(",",ES24.16))') 'C2', J, RHOX(J), &
      (exp(AEXP*RHOX(J)) - exp(AEXP*RHOX(1))) / AEXP, FI(J)
  END DO

  CLOSE(LUOUT)
  WRITE(6, '(A)') 'joshtest: wrote josh_fortran.csv'

CONTAINS

  REAL(8) FUNCTION BNUF(TT)
    REAL(8), INTENT(IN) :: TT
    BNUF = 2.0D0 * HPL * FREQ**3 / CL**2 / (exp(HPL*FREQ/(KB*TT)) - 1.0D0)
  END FUNCTION BNUF

END PROGRAM JOSHTEST
