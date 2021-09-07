PROGRAM test

  implicit none

  integer ihot, ihotstp, imhsf, iths, version, node
  real(8) timehsf
  ihot = 68

  OPEN(IHOT,FILE='fort.68',ACCESS='DIRECT',RECL=8)

  IHOTSTP=1
  READ(IHOT,REC=IHOTSTP) version
  ihotstp = 2
  PRINT *, 'Version = ', version
  READ(IHOT,REC=IHOTSTP) IMHSF
  PRINT *, 'IMHSF is ', IMHSF
  IHOTSTP=3
  READ(IHOT,REC=IHOTSTP) TIMEHSF
  PRINT *, 'TIME is ', TIMEHSF
  IHOTSTP=4
  READ(IHOT, rec=ihotstp) iths
  PRINT *, 'ITHS is ', iths
  IHOTSTP = 5
  read(ihot, rec=ihotstp) node
  print *, 'node = ', node



END PROGRAM
