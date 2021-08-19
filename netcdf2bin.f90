PROGRAM netcdf2bin

  USE netcdf
  IMPLICIT NONE
  include 'netcdf.inc'

  INTEGER i
  INTEGER ncid, status
  INTEGER node, nodeID
  INTEGER recnum

  ! Variables in the hotstart file
  REAL, ALLOCATABLE :: eta1(:), eta2(:), U(:), V(:)
  INTEGER, ALLOCATABLE :: nodecode(:)
  INTEGER IMHS, ITHS, IESTP, NSCOUE, IVSTP, NSCOUV, ICSTP, NSCOUC, IPSTP
  INTEGER IWSTP, NSCOUM, IGEP, NSCOUGE, IGVP, NSCOUGV
  INTEGER IGCP, NSCOUGC, IGPP, IGWP, NSCOUGW
  REAL time


  status = nf90_open(path = "fort.68.nc", mode = nf90_nowrite, ncid = ncid)
  IF (status /= nf90_noerr) THEN
     PRINT *, "Error opening file"
     CALL EXIT(1)
  END if


  ! Read in the integer variables
  ! Commented variables are associated with fort.81 file
  CALL getNetcdfInt(ncid, 'imhs', imhs)
  CALL getNetcdfReal(ncid, 'time', time)
  CALL getNetcdfInt(ncid, 'iths', iths)
  CALL getNetcdfInt(ncid, 'iestp', iestp)
  CALL getNetcdfInt(ncid, 'nscoue', nscoue)
  CALL getNetcdfInt(ncid, 'ivstp', ivstp)
  CALL getNetcdfInt(ncid, 'nscouv', nscouv)
  ! CALL getNetcdfInt(ncid, 'icstp', icstp)
  icstp = 0
  ! CALL getNetcdfInt(ncid, 'nscouc', nscouc)
  nscouc = 0
  CALL getNetcdfInt(ncid, 'ipstp', ipstp)
  CALL getNetcdfInt(ncid, 'iwstp', iwstp)
  CALL getNetcdfInt(ncid, 'nscoum', nscoum)
  CALL getNetcdfInt(ncid, 'igep', igep)
  CALL getNetcdfInt(ncid, 'nscouge', nscouge)
  CALL getNetcdfInt(ncid, 'igvp', igvp)
  CALL getNetcdfInt(ncid, 'nscougv', nscougv)
  ! CALL getNetcdfInt(ncid, 'igcp', igcp)
  igcp = 0
  ! CALL getNetcdfInt(ncid, 'nscougc', nscougc)
  nscougc = 0
  CALL getNetcdfInt(ncid, 'igpp', igpp)
  CALL getNetcdfInt(ncid, 'igwp', igwp)
  CALL getNetcdfInt(ncid, 'nscougw', nscougw)

  ! Get number of nodes and allocate the arrays
  status = nf90_inq_dimid(ncid, 'node', nodeID)
  status = nf90_inquire_dimension(ncid, nodeID, len=node)
  PRINT *, 'Number of nodes = ', node

  ALLOCATE (eta1(node), eta2(node), U(node), V(node), nodecode(node))
  CALL getNetcdfArrayReal(ncid, 'zeta1', node, eta1)
  CALL getNetcdfArrayReal(ncid, 'zeta2', node, eta2)
  CALL getNetcdfArrayReal(ncid, 'u-vel', node, U)
  CALL getNetcdfArrayReal(ncid, 'v-vel', node, V)
  CALL getNetcdfArrayInt(ncid, 'nodecode', node, nodecode)

  !-----------------------------------------------------------------------------------

  ! Write everything to a binary file
  OPEN(unit=1, file='binary', form='unformatted', access='direct', recl=8)

  WRITE(1, rec=1) imhs
  WRITE(1, rec=2) time
  WRITE(1, rec=3) iths

  recnum = 3

  DO i = 1,node
     WRITE(1, rec=recnum+1) eta1(i)
     WRITE(1, rec=recnum+2) eta2(i)
     WRITE(1, rec=recnum+3) U(i)
     WRITE(1, rec=recnum+4) V(i)
     WRITE(1, rec=recnum+5) nodecode(i)

     recnum = recnum + 5
  END DO

  recnum = recnum + 1
  WRITE(1, rec=recnum) IESTP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUE
  recnum = recnum + 1
  WRITE(1, rec=recnum) IVSTP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUV
  recnum = recnum + 1
  WRITE(1, rec=recnum) ICSTP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUC
  recnum = recnum + 1
  WRITE(1, rec=recnum) IPSTP
  recnum = recnum + 1
  WRITE(1, rec=recnum) IWSTP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUM
  recnum = recnum + 1
  WRITE(1, rec=recnum) IGEP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUGE
  recnum = recnum + 1
  WRITE(1, rec=recnum) IGVP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUGV
  recnum = recnum + 1
  WRITE(1, rec=recnum) IGCP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUGC
  recnum = recnum + 1
  WRITE(1, rec=recnum) IGPP
  recnum = recnum + 1
  WRITE(1, rec=recnum) IGWP
  recnum = recnum + 1
  WRITE(1, rec=recnum) NSCOUGW


  CLOSE(1)


END PROGRAM netcdf2bin




SUBROUTINE getNetcdfInt(ncid, name, val)
  USE netcdf
  IMPLICIT NONE
  INTEGER ncid, status, id
  CHARACTER (len = *), INTENT(in) :: name
  INTEGER, INTENT(out) :: val

  status = nf90_inq_varid(ncid, name, id)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error inquiring variable: ", name
     CALL EXIT(1)
  end if

  status = nf90_get_var(ncid, id, val)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error retrieving variable: ", name
     CALL exit(1)
  end if

END SUBROUTINE getNetcdfInt

SUBROUTINE getNetcdfReal(ncid, name, val)
  USE netcdf
  IMPLICIT NONE
  INTEGER ncid, status, id
  CHARACTER (len = *), INTENT(in) :: name
  REAL, INTENT(out) :: val

  status = nf90_inq_varid(ncid, name, id)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error inquiring variable: ", name
     CALL EXIT(1)
  END IF

  status = nf90_get_var(ncid, id, val)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error retrieving variable: ", name
     CALL EXIT(1)
  END if

END SUBROUTINE getNetcdfReal


SUBROUTINE getNetcdfArrayReal(ncid, name, size, val)
  USE netcdf
  IMPLICIT NONE
  INTEGER ncid, status, id, size
  CHARACTER (len = *), INTENT(in) :: name
  REAL, INTENT(out) :: val(size)

  status = nf90_inq_varid(ncid, name, id)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error inquiring array: ", name
     CALL EXIT(1)
  END if
  status = nf90_get_var(ncid, id, val)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error retrieving array: ", name
     CALL EXIT(1)
  END if

END SUBROUTINE getNetcdfArrayReal

SUBROUTINE getNetcdfArrayInt(ncid, name, size, val)
  USE netcdf
  IMPLICIT NONE
  INTEGER ncid, status, id, size
  CHARACTER (len = *), INTENT(in) :: name
  INTEGER, INTENT(out) :: val(size)

  status = nf90_inq_varid(ncid, name, id)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error inquiring array: ", name
     CALL EXIT(1)
  END if
  status = nf90_get_var(ncid, id, val)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error retrieving array: ", name
     CALL EXIT(1)
  END if

END SUBROUTINE getNetcdfArrayInt
