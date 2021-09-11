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