! pftest: standalone probe driver for partition-function interpolation.
!
! Reads probe requests from pf_probes.txt (one per line):
!     H2  <T>
!     BC  <Z> <ION> <T>
!     NEG <Z> <T>
! and writes pf_fortran.csv with the value returned by the ACTUAL
! production code paths (PARTFNH2, U_BC, U_BC_NEG), for cross-checking
! against an independent Python implementation of the same algorithms.
PROGRAM PFTEST
  USE mod_partition_functions, ONLY: U_BC, U_BC_NEG, set_bc_data_dir
  USE mod_atlas_data, ONLY: PARTFNH2, DATADIR
  IMPLICIT NONE

  CHARACTER(len=16)  :: TAG
  CHARACTER(len=256) :: LINE
  REAL(8) :: T, U
  INTEGER :: IZ, ION, IOS, LUNIN, LUNOUT

  DATADIR = '/Users/cconroy/kurucz/atlas12/data/'
  CALL set_bc_data_dir(TRIM(DATADIR))

  OPEN(NEWUNIT=LUNIN,  FILE='pf_probes.txt',  STATUS='OLD', ACTION='READ')
  OPEN(NEWUNIT=LUNOUT, FILE='pf_fortran.csv', STATUS='REPLACE', ACTION='WRITE')

  DO
    READ(LUNIN, '(A)', IOSTAT=IOS) LINE
    IF (IOS .NE. 0) EXIT
    IF (LEN_TRIM(LINE) .EQ. 0) CYCLE
    READ(LINE, *) TAG
    SELECT CASE (TRIM(TAG))
    CASE ('H2')
      READ(LINE, *) TAG, T
      U = PARTFNH2(T)
      WRITE(LUNOUT, '(A,",",I3,",",I2,",",ES23.15,",",ES23.15)') 'H2', 0, 0, T, U
    CASE ('BC')
      READ(LINE, *) TAG, IZ, ION, T
      U = U_BC(IZ, ION, T)
      WRITE(LUNOUT, '(A,",",I3,",",I2,",",ES23.15,",",ES23.15)') 'BC', IZ, ION, T, U
    CASE ('NEG')
      READ(LINE, *) TAG, IZ, T
      U = U_BC_NEG(IZ, T)
      WRITE(LUNOUT, '(A,",",I3,",",I2,",",ES23.15,",",ES23.15)') 'NEG', IZ, 0, T, U
    END SELECT
  END DO
  CLOSE(LUNIN)
  CLOSE(LUNOUT)
END PROGRAM PFTEST
