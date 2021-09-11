program dumpoutput
  use netcdf
  implicit none

  real(8), allocatable :: eta2(:)
  INTEGER i
  INTEGER ncid, status
  INTEGER node, nodeID, nele, neleID
  INTEGER recnum
  INTEGER nodecodeID
  integer eta2id

  status = nf90_open(path='fort.68.nc', mode = nf90_nowrite, ncid = ncid)
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

  ALLOCATE (eta2(node))
  status = nf90_inq_varid(ncid, 'zeta2', eta2id)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error inquiring array: zeta2"
     CALL EXIT(1)
  END if
  status = nf90_get_var(ncid, eta2id, eta2)
  if (status /= nf90_noerr) THEN
     PRINT *, "Error retrieving array: zeta2"
     CALL EXIT(1)
  END if

  ! write to text file
  open(63, file='fort.63.dump')

2453              FORMAT(2X,I8,2X,E16.8E3)
  do i = 1,node
     write(63, 2453) i, eta2(i)
  end do

end program dumpoutput
