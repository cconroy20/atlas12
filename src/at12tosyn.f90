PROGRAM at12tosyn

  !routine to convert an ATLAS12 output file into a SYNTHE input file
  !calling sequence: at12tosyn.exe infile outfile
  
  INTEGER :: i,nlines
  CHARACTER(150) :: ifile,ofile
  CHARACTER(99) :: str
  CHARACTER(2)  :: str2

  CALL GETARG(1,ifile)
  CALL GETARG(2,ofile)

  !open ATLAS12 output file
  OPEN(33,file=TRIM(ifile),ACTION='READ',STATUS='OLD')

  !open SYNTHE input file
  OPEN(34,FILE=TRIM(ofile),ACTION='WRITE',STATUS='REPLACE')
  !add the SYNTHE cards
  WRITE(34,'("SURFACE FLUX")') 
  WRITE(34,'("ITERATIONS 1 PRINT 2 PUNCH 2")') 
  WRITE(34,'("CORRECTION OFF")') 
  WRITE(34,'("PRESSURE OFF")') 
  WRITE(34,'("READ MOLECULES")') 
  WRITE(34,'("MOLECULES ON")') 

  !re-write the ATLAS12 header
  DO i=1,22
     READ(33,'(A99)') str
     WRITE(34,'(A99)') str
  ENDDO

  !burn the extra ATLAS12 stuff
  DO i=1,22
     READ(33,'(A99)') str
  ENDDO

  !write the model structure
  READ(33,'(A99)') str
  WRITE(34,'(A99)') str

  READ(str(12:14),'(I2)') nlines

  DO i=1,nlines+2
     READ(33,'(A99)') str
     WRITE(34,'(A99)') str
  ENDDO

  CLOSE(33)
  CLOSE(34)

END PROGRAM at12tosyn
