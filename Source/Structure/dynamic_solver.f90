subroutine dynamic_solver(the_model_structure, settings, boolean_output_write)
! Importing classes
  use my_math_structure, only: cross, eye, outer
  use my_constants_structure, only: e1, e2, e3
  use class_model_structure, only: model_structure, step
  use solver_variables
  use mkl_dss
  use class_sparse_matrix
  use solver_functions
  use my_FileIO_functions
  
  implicit none
  
  real(kind=8), allocatable :: fqext0(:)
  real(kind = 8), allocatable :: statOut(:)

  type(model_structure), intent(inout) :: the_model_structure
  type(step), intent(inout) :: settings
  type(MKL_DSS_HANDLE)	:: handle
  
  logical :: boolean_sparse = .TRUE., boolean_external_force = .FALSE.
  logical, intent(in) :: boolean_output_write

  real(kind = 8) :: max_de
  
  ! Variables to test a modified time step
  integer :: iterations !iterations needed to converge
  
!< Switching warnings off. PARDISO compatibility problem with new version OPENMP! Maybe resolved in newer version
  !call kmp_set_warnings_off() 
  
  convert_job(1) = 1
  convert_job(2) = 1
  convert_job(3) = 1
  convert_job(4) = 2

! Dimensioning result vectors for dynamic calculation
  nMatSize = 2*tncoordinates + trconstraints
    
  allocate(deltaxvector(nMatSize))
  allocate(rhs_unbalance_force(nMatSize))
  allocate(fqext0(tncoordinates))
  
! create direct sparse solver matrix
  error = dss_create(handle, MKL_DSS_DEFAULTS)
  error = dss_define_structure(handle, MKL_DSS_NON_SYMMETRIC, sparse_smatrix%csr_rowIndices, nMatSize, nMatSize, sparse_smatrix%csr_columns, size(sparse_smatrix%csr_columns) );
  error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
  
  if (error /= 0) then
    boolean_sparse = .FALSE.
  endif

! Initialization of field variables
  q0        = the_model_structure%q_t
  v0        = the_model_structure%qdot_t
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
  
  print '(2X,A37,1X,F10.5)', "Starting dynamic calculation at time:", the_model_structure%time_t
  
! Loop over time steps
  !do while (istep .le. inos)
  do while (settings%totalt .gt. time)
      
!!------------------------------------------------------------------------------------------------------
!! Adaptation of the time step based on the iterations of the previous time step (aim for 3 iterations)
!!------------------------------------------------------------------------------------------------------
!if (settings%optional_adaptive_deltat == 1) then      
!  if (time .ge. settings%deltat) then
!    if (iterations .ge. 6) then
!        the_model_structure%deltat = max((0.7d0)*settings%deltat ,1.0d-6)
!        settings%deltat            = max((0.7d0)*settings%deltat ,1.0d-6)
!        print*, 'delta t adapted to ', settings%deltat 
!    elseif(iterations .lt. 5) then
!        the_model_structure%deltat = min(1.3d0*settings%deltat ,3.0d1)
!        settings%deltat            = min(1.3d0*settings%deltat ,3.0d1)
!        print*, 'delta t adapted to ', settings%deltat 
!    end if       
!  end if  
!end if 
!!------------------------------------------------------------------------------------------------------
!      
! update time step
    istep = istep + 1
    time  = time + settings%deltat
    
! Get model loads
    if (time .eq. 0.0d0) then
      call the_model_structure%modelloads(time, ml_t, sl_t)
    else
      call the_model_structure%modelloads(time-settings%deltat/2, ml_t, sl_t)
    end if
     
! Prescribed external force vector: fqext0(t) = fqext0(n) + xhi*(fqext0(n+1) - fqext0(n))
    if (boolean_external_force) then
      fqext0  = the_model_structure%fqext0_t + (the_model_structure%fqext0 - the_model_structure%fqext0_t) / (settings%totalt - the_model_structure%time_t) * (time - the_model_structure%time_t)
    end if
  
! Get model boundaries at current time increment
    call the_model_structure%modelboundaries(time)
     
! Configuration for positions and velocities at current time step
    q1 = q0
    v1 = v0
     
! Start Newton iteration
    iteration = 0
    residuum  = 1.0d0
    do while((residuum > settings%tolerance).and.(iteration < settings%iterationlimit))
      iteration = iteration+1

! computing f and k for the whole model.
      call sparse_smatrix%startAssembly()  
      call the_model_structure%modelfk(q0, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q1, v0, v1)
        
! computing g and H = dg for the whole model at q0 and q1.
      call the_model_structure%modelgdg(q1, q0, sparse_smatrix, indicesg, indicesq, lambda)
        
! computing constraint forces corresponding to q
      call the_model_structure%modelfL(lambda,sparse_smatrix)
        
! construction of the residual vector
      rhs_unbalance_force(indicesq) = the_model_structure%fqint + the_model_structure%fqdyn + the_model_structure%fqL - the_model_structure%fqgrav - the_model_structure%fqext - fqext0
      rhs_unbalance_force(indicesv) = the_model_structure%fv
      rhs_unbalance_force(indicesg) = the_model_structure%g
      
! Solving linear system of equation<b
#if !B05_ADDON
      call solver(nMatSize, handle, sparse_smatrix%csr_values, sparse_smatrix%csr_columns, sparse_smatrix%csr_rowIndices, rhs_unbalance_force, 1, deltaxvector, boolean_sparse, statOut, .FALSE.)
#endif      
      residuum = dsqrt(dot_product(deltaxvector, deltaxvector))/dsqrt(dot_product(q1, q1)+dot_product(v1, v1)+dot_product(lambda, lambda))
      
      !print*, dsqrt(dot_product(rhs_unbalance_force, rhs_unbalance_force)), dsqrt(dot_product(deltaxvector, deltaxvector))
      
      if (boolean_output_write) then
        print '(5X,A9,1X,E10.2E2,1X,A17,1X,E10.2E3)', 'residuum:', residuum, 'convergence work:',abs(dot_product(rhs_unbalance_force, deltaxvector))
      end if
      
      if (isnan(residuum)) then
        the_model_structure%boolean_abort = .True.
        exit
      endif
      
! additive actualization of generalized displacements and Lagrange multiplier
      q1     = q1     + deltaxvector(indicesq)
      v1     = v1     + deltaxvector(indicesv)
      lambda = lambda + deltaxvector(indicesg)

    end do
        
! Abort criteria
    if (the_model_structure%boolean_abort) exit
    if (residuum > settings%tolerance) then
      print*, 'solution failed to converge at time', time
      the_model_structure%boolean_abort = .True.
      exit
    end if

! calculate invariants of time current step
    call the_model_structure%modelinvariants(q1, v1)

! Actualization of the old soltion     
    q0    = q1
    v0    = v1

! writing in result file
    if (boolean_output_write) then
      call the_model_structure%output_write(settings%step,time, time - the_model_structure%time_t, istep, settings%i_simutype, iteration, q1, v1, lambda, the_model_structure%stress)
    end if
     
    print '(5X,A5,1X,F10.5,1X,A11,1X,I2,1X,A9,1X,E10.2E3)','time:',time,'iterations:',iteration,'residuum:',residuum

! iterations need for time step 
    iterations = iteration

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
  the_model_structure%fqext0_t       = fqext0
  the_model_structure%fqint_t        = the_model_structure%fqint
  the_model_structure%fqdyn_t        = the_model_structure%fqdyn
  the_model_structure%fqL_t          = the_model_structure%fqL
  the_model_structure%fv_t           = the_model_structure%fv

  if (settings%output_flag == 1) then
    print *, "Output: Writing matrices..."
    
    call writeIntVectorToFile(indicesq, 9999,'indicesq.dres')
    call writeIntVectorToFile(indicesg, 9999,'indicesg.dres')
    call writeIntVectorToFile(indicesv, 9999,'indicesv.dres')
      
    call writeIntVectorToFile(sparse_smatrix%csr_columns,   9999,'csr_columns_Kaug.dres')
    call writeIntVectorToFile(sparse_smatrix%csr_rowIndices,9999,'csr_rowIndices_Kaug.dres')
    call writeRealVectorToFile(sparse_smatrix%csr_values,   9999,'csr_values_Kaug.dres')
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
  deallocate(rhs_unbalance_force)
  deallocate(fqext0)
  return
end subroutine dynamic_solver
