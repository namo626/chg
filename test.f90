PROGRAM test

  implicit none

  integer ihot, ihotstp, imhsf, iths
  real timehsf
  ihot = 68

  OPEN(IHOT,FILE='fort.68',ACCESS='DIRECT',RECL=8)

  IHOTSTP=1
  READ(IHOT,REC=IHOTSTP) IMHSF
  PRINT *, 'IMHSF is ', IMHSF
  IHOTSTP=2
  READ(IHOT,REC=IHOTSTP) TIMEHSF
  PRINT *, 'TIME is ', TIMEHSF
  IHOTSTP=3
  READ(IHOT, rec=ihotstp) iths
  PRINT *, 'ITHS is ', iths



END PROGRAM
