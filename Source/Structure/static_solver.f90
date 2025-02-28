subroutine static_solver(the_model_structure, settings, boolean_output_write)
! Importing classes
  use my_math_structure, only: cross, eye
  use class_model_structure, only: model_structure, step
  use solver_variables
  use mkl_dss
  use class_sparse_matrix
  use solver_functions
  use my_FileIO_functions
  
  implicit none
  type(model_structure), intent(inout) :: the_model_structure
  type(step), intent(in)     :: settings
  type(MKL_DSS_HANDLE)	:: handle
  
  logical :: boolean_sparse = .TRUE., boolean_external_force = .FALSE.
  logical :: boolean_determinant = .TRUE.
  logical, intent(in) :: boolean_output_write

  real(kind = 8), allocatable :: statOut(:)
  real(kind = 8), dimension(:), allocatable :: csr_values_stat, csr_values_kqq 
  real(kind=8), allocatable :: fqgrav0(:), fqext0(:)
  
  integer, dimension(:), allocatable :: csr_columns_stat, csr_rowIndices_stat, indices_qg, csr_compr_source_stat, csr_columns_kqq , csr_rowIndices_kqq , csr_compr_source_kqq
  
  allocate(statOut(5))
  statOut = 0.0d0

!< Switching warnings off. PARDISO compatibility problem with new version OPENMP! Maybe resolved in newer version
  !call kmp_set_warnings_off() 
  
  convert_job(1) = 1
  convert_job(2) = 1
  convert_job(3) = 1
  convert_job(4) = 2

! Dimensioning result vectors for static calculation
  nMatSize = 2*tncoordinates + trconstraints
  
! allocating solution vector  
  allocate(deltaxvector(nMatSize-tncoordinates))
  allocate(rhs_unbalance_force(nMatSize-tncoordinates))
  allocate(fqgrav0(tncoordinates))
  allocate(fqext0(tncoordinates))
  allocate(indices_qg(nMatSize-tncoordinates))
  indices_qg = reshape((/indicesq,indicesg/),(/nMatSize-tncoordinates/))

! Creating and extracting sparse matrix from dynamic system matrix for static analysis
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, csr_values_stat, csr_columns_stat, csr_rowIndices_stat, csr_compr_source_stat, indices_qg, indices_qg)
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, csr_values_kqq , csr_columns_kqq , csr_rowIndices_kqq , csr_compr_source_kqq , indicesq,   indicesq)  
  
! Getting new solver handle for static iteration matrix
  !print *, "Reordering the static iteration matrix..."
  error = dss_create(handle, MKL_DSS_DEFAULTS)
  error = dss_define_structure(handle, MKL_DSS_NON_SYMMETRIC, csr_rowIndices_stat, nMatSize-tncoordinates, nMatSize-tncoordinates, csr_columns_stat, size(csr_columns_stat))
  error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
      
  if (error /= 0) then
    boolean_sparse = .FALSE.
  endif

! Initialization of field variables
  q1        = the_model_structure%q_t
  v1        = the_model_structure%qdot_t 
  lambda(:) = the_model_structure%lambda_t
  
  the_model_structure%deltat = settings%deltat
  time  = the_model_structure%time_t
  inos  = IDNINT((settings%totalt - the_model_structure%time_t)/the_model_structure%deltat) - 1
  istep = 0
  
  fqext0(:) = the_model_structure%fqext0_t
  if (norm2(the_model_structure%fqext0) .gt. 0.0d0) boolean_external_force = .TRUE.
  
  the_model_structure%fqgrav(:) = 0.0d0
  if (settings%gravityflag == 1) then
    call the_model_structure%modelfgrav()
  end if
  fqgrav0 = the_model_structure%fqgrav
  
  print '(2X,A36,1X,F10.5)', "Starting static calculation at step:", the_model_structure%time_t
  
! Loop over time steps
  do while (istep .le. inos)

! update time step
    istep = istep + 1
    time  = time + the_model_structure%deltat
       
! Get model loads
    call the_model_structure%modelloads(time, ml_t, sl_t)
    ! gravity forces
    if (settings%gravityflag == 1) then
      the_model_structure%fqgrav = fqgrav0*(time - the_model_structure%time_t) / (settings%totalt - the_model_structure%time_t)
    end if
    ! Prescribed external force vector: fqext0(t) = fqext0(n) + xhi*(fqext0(n+1) - fqext0(n))
    if (boolean_external_force) then
      fqext0  = the_model_structure%fqext0_t + (the_model_structure%fqext0 - the_model_structure%fqext0_t) / (settings%totalt - the_model_structure%time_t) * (time - the_model_structure%time_t)
    end if
    
! Get model boundaries at current time increment
    call the_model_structure%modelboundaries(time)
     
! Start Newton iteration
    iteration = 0
    residuum  = 1.0d0
    q0 = q1
    do while((residuum > settings%tolerance).and.(iteration < settings%iterationlimit))
      iteration = iteration+1
		
! get problem matrix and residual vector
      call sparse_smatrix%startAssembly()
      call the_model_structure%modelfk(q1, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q1)
		
! computing g and dg for the whole model at q0 and q1.
      call the_model_structure%modelgdg(q1, q1, sparse_smatrix, indicesg, indicesq, lambda, q0)
	
! computing constraint forces corresponding to q
      call the_model_structure%modelfL(lambda,sparse_smatrix)
        
! construction of the residual vector
      rhs_unbalance_force(indicesq)      = the_model_structure%fqint + the_model_structure%fqdyn + the_model_structure%fqL - the_model_structure%fqext - the_model_structure%fqgrav - fqext0
      rhs_unbalance_force(indicesg_stat) = the_model_structure%g

      if (boolean_sparse) then
        sparse_smatrix%csr_values(csr_compr_source_kqq) = sparse_smatrix%csr_values(csr_compr_source_kqq)*2.0d0
        csr_values_stat = 0.0d0
        csr_values_stat = sparse_smatrix%csr_values(csr_compr_source_stat)
      end if    

! Solving linear system of equations
#if !B05_ADDON
      call solver(nMatSize-tncoordinates, handle, csr_values_stat, csr_columns_stat, csr_rowIndices_stat, rhs_unbalance_force, 1, deltaxvector, boolean_sparse, statOut, boolean_determinant)
#endif     
      residuum = dsqrt(dot_product(deltaxvector, deltaxvector))/dsqrt(dot_product(q1, q1)+dot_product(v1, v1)+dot_product(lambda, lambda))
      if (boolean_output_write .eqv. .TRUE.) then
        print '(5X,A9,1X,E10.2E2,1X,A18,1X,E10.2E3)', 'residuum:', residuum, 'unbalanced forces:',dot_product(rhs_unbalance_force, deltaxvector)
      end if
        
      if (isnan(residuum)) then
        the_model_structure%boolean_abort = .True.
        exit
      endif
      
! additive actualization of generalized displacements and Lagrange multiplier
      q1     = q1  + deltaxvector(indicesq)
      v1     = v1
      lambda = lambda + deltaxvector(indicesg_stat)

    end do

! abort criteria
    if (the_model_structure%boolean_abort) exit
    if (residuum > settings%tolerance) then
      if (settings%iterationlimit /= 1) then
        print*, 'solution failed to converge at step ', iStep, 'time', time
        the_model_structure%boolean_abort = .True.
      end if 
      exit
    end if

! calculate invariants of time current step
    call the_model_structure%modelinvariants(q1, v1)

! writing in result file
    if (boolean_output_write) then
      call the_model_structure%output_write(settings%step,time, time - the_model_structure%time_t, istep, settings%i_simutype, iteration, q1, v1, lambda, the_model_structure%stress)
    end if
    
    print '(5X,A5,1X,F10.3,1X,A11,1X,I2,1X,A9,1X,E10.2E3)','time:',time,'iterations:',iteration,'residuum:',residuum
    
  end do

  if (the_model_structure%boolean_abort) then
    print '(5X,A23,1X,F10.5,1X,A11,1X,I2,1X,A9,1X,E10.2E3)','no convergence at time:',time,'iterations:',iteration,'residuum:',residuum
    goto 999
  end if
  
! Updating model variables
  the_model_structure%q_t            = q1
  the_model_structure%qdot_t         = v1
  the_model_structure%lambda_t       = lambda
  the_model_structure%time_t         = time
  
  the_model_structure%kenergy_t      = the_model_structure%kenergy
  the_model_structure%strainenergy_t = the_model_structure%strainenergy
  the_model_structure%gravenergy_t   = the_model_structure%gravenergy
  the_model_structure%penergy_t      = the_model_structure%penergy
  the_model_structure%tenergy_t      = the_model_structure%tenergy
  
  the_model_structure%lmomentum_t    = the_model_structure%lmomentum
  the_model_structure%amomentum_t    = the_model_structure%amomentum
  
  the_model_structure%fqext_t        = the_model_structure%fqext
  the_model_structure%fqint_t        = the_model_structure%fqint
  the_model_structure%fqdyn_t        = 0.0d0
  the_model_structure%fqL_t          = the_model_structure%fqL
  the_model_structure%fv_t           = 0.0d0
  the_model_structure%fqext0_t       = fqext0
  
  if (settings%output_flag == 1) then
    print *, "Output: Writing matrices..."

    call writeIntVectorToFile(indicesq, 9999,'indicesq.dres')
    call writeIntVectorToFile(indicesg_stat, 9999,'indicesg.dres')
      
    call writeIntVectorToFile(csr_columns_stat,   9999,'csr_columns_Kstat.dres')
    call writeIntVectorToFile(csr_rowIndices_stat,9999,'csr_rowIndices_Kstat.dres')
    call writeRealVectorToFile(csr_values_stat,   9999,'csr_values_Kstat.dres')
  end if   
  
999 continue
    
  allocate(temp_check(1))
  if (the_model_structure%boolean_abort) then
    temp_check = 0.0d0
    call writeRealVectorToFile(temp_check,9999,'check.log')
  else
    temp_check = dsqrt(dot_product(the_model_structure%q_t,the_model_structure%q_t))
    call writeRealVectorToFile(temp_check,9999,'check.log')
  end if
  deallocate(temp_check)
  
  deallocate(deltaxvector)
  deallocate(indices_qg)
  deallocate(rhs_unbalance_force)
  deallocate(csr_rowIndices_stat)
  deallocate(csr_columns_stat)
  deallocate(csr_values_stat)
  deallocate(csr_compr_source_stat)
  deallocate(fqext0)
  deallocate(fqgrav0)
  
  return
end subroutine static_solver
