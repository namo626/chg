PROGRAM netcdf2bin
  ! This program converts a hotstart file in NetCDF format to a binary format
  ! that can be read by DG-SWEM. Some variables were left out or manually
  ! defined because this program was written for a specific ADCIRC hotstart file;
  ! you can manually call them as necessary (e.g. by looking at ncdump) in the
  ! appropriate record order.
  ! This program was tested with ifort and netcdf 4.6.2. There seems to be a problem
  ! with a newer version of netcdf.

  ! COMPILE
  ! ifort netcdf2bin.f90 -o netcdf2bin -I{NETCDF_INC} -L{NETCDF_LIB} -lnetcdf -lnetcdff

  ! USAGE
  ! netcdf2bin INPUT OUTPUT
  ! Example: netcdf2bin fort.68.nc fort.68

  USE netcdf
  IMPLICIT NONE

  CHARACTER (100) :: input, output


  INTEGER i
  INTEGER ncid, status
  INTEGER node, nodeID, nele, neleID
  INTEGER recnum
  INTEGER nodecodeID

  ! Variables in the hotstart file
  REAL(8), ALLOCATABLE :: eta1(:), eta2(:), U(:), V(:)
  INTEGER, ALLOCATABLE :: nodecode(:), noff(:)
  INTEGER IMHS, ITHS, IESTP, NSCOUE, IVSTP, NSCOUV, ICSTP, NSCOUC, IPSTP
  INTEGER IWSTP, NSCOUM, IGEP, NSCOUGE, IGVP, NSCOUGV
  INTEGER IGCP, NSCOUGC, IGPP, IGWP, NSCOUGW
  REAL(8) time

  IF (command_argument_count().ne.2) THEN
     PRINT *, 'Error: need 2 arguments'
     CALL EXIT(1)
  END IF

  CALL get_command_argument(1, input)
  CALL get_command_argument(2, output)

  status = nf90_open(path = input, mode = nf90_nowrite, ncid = ncid)
  IF (status /= nf90_noerr) THEN
     PRINT *, "Error opening hotstart file"
     CALL EXIT(1)
  END if

  ! Get number of nodes and allocate the arrays
  status = nf90_inq_dimid(ncid, 'node', nodeID)
  status = nf90_inquire_dimension(ncid, nodeID, len=node)
  PRINT *, 'Number of nodes = ', node

  status = nf90_inq_dimid(ncid, 'nele', neleID)
  status = nf90_inquire_dimension(ncid, neleID, len=nele)
  PRINT *, 'Number of elements = ', nele

  ALLOCATE (eta1(node), eta2(node), U(node), V(node))
  ALLOCATE (nodecode(node), noff(nele))


  PRINT *, 'Start retrieving variables...'
  ! Add required variables as necessary

  ! IM 111112 not in dgswem, so set to zero
  ! CALL getNetcdfInt(ncid, 'imhs', imhs)
  imhs = 0
  CALL getNetcdfReal(ncid, 'time', time)
  PRINT *, 'Time = ', time
  CALL getNetcdfInt(ncid, 'iths', iths)
  PRINT *, 'ITHS = ', iths
  CALL getNetcdfInt(ncid, 'iestp', iestp)
  CALL getNetcdfInt(ncid, 'nscoue', nscoue)
  PRINT *, 'NSCOUE = ', nscoue
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
  PRINT *, 'NSCOUGW = ', nscougw


  CALL getNetcdfArrayReal(ncid, 'zeta1', node, eta1)
  call getNetcdfArrayReal(ncid, 'zeta2', node, eta2)
  CALL getNetcdfArrayReal(ncid, 'u-vel', node, U)
  CALL getNetcdfArrayReal(ncid, 'v-vel', node, V)
  CALL getNetcdfArrayInt(ncid, 'nodecode', node, nodecode)
  CALL getNetcdfArrayInt(ncid, 'noff', nele, noff)

  PRINT *, 'Finished retrieving netCDF variables.'

  !-----------------------------------------------------------------------------------

  ! Write everything to a binary file
  OPEN(unit=1, file=output,  access='direct', recl=8)
  PRINT *, 'Writing results to binary...'
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
     WRITE(1, rec=recnum+6) noff(i)

     recnum = recnum + 6
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

  INCLUDE 'netcdf2bin.inc'

END SUBROUTINE getNetcdfInt

SUBROUTINE getNetcdfReal(ncid, name, val)
  USE netcdf
  IMPLICIT NONE
  INTEGER ncid, status, id
  CHARACTER (len = *), INTENT(in) :: name
  REAL(8), INTENT(out) :: val

  INCLUDE 'netcdf2bin.inc'

END SUBROUTINE getNetcdfReal


SUBROUTINE getNetcdfArrayReal(ncid, name, size, val)
  USE netcdf
  IMPLICIT NONE
  INTEGER ncid, status, id, size
  CHARACTER (len = *), INTENT(in) :: name
  REAL(8), INTENT(out) :: val(size)

  INCLUDE 'netcdf2bin.inc'

END SUBROUTINE getNetcdfArrayReal

SUBROUTINE getNetcdfArrayInt(ncid, name, size, val)
  USE netcdf
  IMPLICIT NONE
  INTEGER ncid, status, id, size
  CHARACTER (len = *), INTENT(in) :: name
  INTEGER, INTENT(out) :: val(size)

  INCLUDE 'netcdf2bin.inc'

END SUBROUTINE getNetcdfArrayInt
