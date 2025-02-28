subroutine invariants(the_model_structure, settings, boolean_output_write)
  use class_model_structure, only: model_structure, step
  use solver_variables
  use class_sparse_matrix
  use solver_functions
  use my_FileIO_functions
  implicit none
  
  type(model_structure), intent(inout) :: the_model_structure
  type(step), intent(in) :: settings
  logical, intent(in) :: boolean_output_write
  integer :: io_error, nsteps
  
!< Initialization of field variables
  ml_t   = 0.0d0
  sl_t   = 0.0d0
  q1     = 0.0d0
  v1     = 0.0d0
  lambda = 0.0d0
  
  the_model_structure%fqgrav(:) = 0.0d0
  if (settings%gravityflag == 1) then
    call the_model_structure%modelfgrav()
  end if

  open(unit = 9999, file = the_model_structure%output_filename // '_q.dres', status = 'old', iostat = io_error)
  nsteps = 0
  read(9999, *, iostat = io_error)
  do 
    if (io_error .ne. 0) exit
    read(9999, *, iostat = io_error)
    nsteps = nsteps + 1
  end do
  close(unit = 9999)
  nsteps = nsteps - 1
  
  !< Loop over all steps
  open(unit = 9998, file = the_model_structure%output_filename // '_q.dres', status = 'old', iostat = io_error)
  if (io_error .ne. 0) then
      print*, 'Warning: error in reading coordinates file...'
      call exit(-1)
  end if
  
  open(unit = 9999, file = the_model_structure%output_filename // '_v.dres', status = 'old', iostat = io_error)
  if (io_error .ne. 0) then
      print*, 'Warning: no velocity file found...'
  end if
  
  istep = 0
  do while (istep .le. nsteps)
    istep = istep + 1
    
    !< read kinematic from solution file
    read(9998, *, iostat = io_error) q1
    if (io_error .ne. 0) then
      the_model_structure%boolean_abort = .TRUE.
      print*, 'Warning: error in reading coordinates at step', istep
      exit
    end if
    
    read(9999, *, iostat = io_error) v1
    if (io_error .ne. 0) then
      v1 = 0.0d0
    end if
    
    !< calculate potential energy of elements
    call sparse_smatrix%startAssembly()
    call the_model_structure%modelfk(q1, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q1, v1, v1)
    
    !< calculate invariants of time current step
    call the_model_structure%modelinvariants(q1, v1)

    !< writing in result file
    if (boolean_output_write .eqv. .TRUE.) then
      call the_model_structure%output_write(istep, 0.0d0, 0.0d0, istep, 0, 0, q1, v1, the_model_structure%lambda_t, the_model_structure%stress, 6)
      call the_model_structure%output_write(istep, 0.0d0, 0.0d0, istep, 0, 0, q1, v1, the_model_structure%lambda_t, the_model_structure%stress, 7)
      call the_model_structure%output_write(istep, 0.0d0, 0.0d0, istep, 0, 0, q1, v1, the_model_structure%lambda_t, the_model_structure%stress, 8)
    end if
     
    print '(5X,A13,1X,I4)','reading step:', istep

  end do
  
  close(unit = 9998)
  close(unit = 9999)
  
  if (settings%output_flag == 1) then
    print *, "Output: Writing matrices..."
    
    call writeIntVectorToFile(indicesq, 9999,'indicesq.dres')
    call writeIntVectorToFile(indicesg, 9999,'indicesg.dres')
    call writeIntVectorToFile(indicesv, 9999,'indicesv.dres')
      
    call writeIntVectorToFile(sparse_smatrix%csr_columns,   9999,'csr_columns_Kaug.dres')
    call writeIntVectorToFile(sparse_smatrix%csr_rowIndices,9999,'csr_rowIndices_Kaug.dres')
    call writeRealVectorToFile(sparse_smatrix%csr_values,   9999,'csr_values_Kaug.dres')
  end if  

  allocate(temp_check(1))
  if (the_model_structure%boolean_abort) then
    temp_check = 0.0d0
    call writeRealVectorToFile(temp_check,9999,'check.log')
  else
    temp_check = dsqrt(dot_product(q1,q1))
    call writeRealVectorToFile(temp_check,9999,'check.log')
  end if
  deallocate(temp_check)
  
  return
end subroutine invariants
