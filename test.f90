PROGRAM test

  implicit none

  integer ihot, ihotstp, imhsf, timehsf
  ihot = 68

  OPEN(IHOT,FILE='fort.68',ACCESS='DIRECT',RECL=8)

  IHOTSTP=1
  READ(IHOT,REC=IHOTSTP) IMHSF
  IHOTSTP=2
  READ(IHOT,REC=IHOTSTP) TIMEHSF
  ! debug
  PRINT *, 'TIME is ', TIMEHSF
  IHOTSTP=3

      
END PROGRAM
