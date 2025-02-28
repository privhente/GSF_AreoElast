subroutine static_solver_arc(the_model_structure, settings, boolean_output_write)
! Declaring variables
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
  TYPE(MKL_DSS_HANDLE)	:: handle
  logical :: boolean_sparse = .TRUE.
  logical :: boolean_determinant = .FALSE.
  logical, intent(in) :: boolean_output_write
  
  real(kind = 8), dimension(:), allocatable :: sparse_csr_values_K_stat, sparse_csr_values_kqq 
  integer, dimension(:), allocatable :: sparse_csr_columns_K_stat, sparse_csr_rowIndices_K_stat, indices_qg, sparse_csr_compr_source_K_stat, sparse_csr_columns_kqq , sparse_csr_rowIndices_kqq , sparse_csr_compr_source_kqq

! solution variables  
  real(kind = 8), allocatable :: delta_q(:), lambdan(:), delta_lambda(:)

! variables for the arc-length method
  real(kind = 8), allocatable :: delta_x_solver(:,:), rhs_solver(:,:), delta_q_r(:), delta_q_g(:), qn(:), D_q(:), D_qn(:) 
  real(kind = 8), allocatable :: delta_lambda_g(:), delta_lambda_r(:)
  real(kind = 8) :: ffun, alfa, alfan, delta_alfa, D_alfa, f_alfa, dL_arc, dL_arc0, dL_arc_n, DL_arc_min
  real(kind = 8), allocatable :: f_q(:)
  real(kind = 8) :: sign, stiffness_parameter, I_n
  real(kind = 8) :: a1, a2, a3, a4, a5, delta_alfa_1, delta_alfa_2, cos1, cos2, crit_u
  
  integer, parameter :: bufLen = 20
  character*15 statIn
  integer buff(bufLen)
  real(kind = 8), allocatable :: statOut(:)

!< Switching warnings off. PARDISO compatibility problem with new version OPENMP! Maybe resolved in newer version
  !call kmp_set_warnings_off() 
  
  convert_job(1) = 1
  convert_job(2) = 1
  convert_job(3) = 1
  convert_job(4) = 2
  
! Dimensioning result vectors for static calculation
  nq       = tncoordinates
  nc       = trconstraints
  nMatSize = nq + nc
  
! allocating solution vector  
  allocate(statOut(5))
  statOut = 0.0d0
    
  allocate(delta_lambda(nc))
  allocate(delta_lambda_r(nc))
  allocate(delta_lambda_g(nc))
  allocate(lambdan(nc))
  
  allocate(delta_q(nq))
  allocate(delta_q_r(nq))
  allocate(delta_q_g(nq))
  
  allocate(delta_x_solver(nMatSize,2))
  allocate(deltaxvector(nMatSize))
  
  allocate(rhs_unbalance_force(nMatSize))
  allocate(rhs_solver(nMatSize,2))
  
  allocate(qn(nq))
  allocate(D_q(nq))
  allocate(D_qn(nq))
  allocate(f_q(nq))
  
  allocate(indices_qg(nMatSize))
  indices_qg = reshape((/indicesq,indicesg/),(/nMatSize/))

! Creating and extracting sparse matrix from dynamic system matrix for static analysis
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, sparse_csr_values_K_stat, sparse_csr_columns_K_stat, sparse_csr_rowIndices_K_stat, sparse_csr_compr_source_K_stat, indices_qg, indices_qg)
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, sparse_csr_values_kqq , sparse_csr_columns_kqq , sparse_csr_rowIndices_kqq , sparse_csr_compr_source_kqq , indicesq, indicesq)  
  
! Getting new solver handle for static iteration matrix
  error = dss_create(handle, MKL_DSS_DEFAULTS)
  error = dss_define_structure(handle, MKL_DSS_NON_SYMMETRIC, sparse_csr_rowIndices_K_stat, nMatSize, nMatSize, sparse_csr_columns_K_stat, size(sparse_csr_columns_K_stat))
  error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
     
  if (error /= 0) then
    boolean_sparse = .FALSE.
  endif

! Initialization of field variables
  q0        = the_model_structure%q_t
  v0        = the_model_structure%qdot_t 
  lambda(:) = the_model_structure%lambda_t
  lambdan(:)= lambda(:)
  inos      = settings%totalt/settings%deltat

  alfan     = 0.0d0
  qn        = q0
  D_q       = 0.0d0
  D_qn      = 0.0d0
  DL_arc0   = settings%deltat
  
  DL_arc_min= DL_arc0
  DL_arc_n  = DL_arc0
  dL_arc    = DL_arc0
  I_n       =  settings%number_desired_iteration
  iteration = I_n
  
  print '(2X,A43)', "Starting static load-control calculation..."
  
! Loop over time steps
  time   = 1.0d0
  do while( alfan <= settings%totalt ) 
    inos = inos + 1

! Update unknowns
    q1     = qn
    D_q    = D_qn
    lambda = lambdan
    
! Get model_structure loads
     call the_model_structure%modelloads(time, ml_t, sl_t)

! Get model_structure boundaries at current time increment
     call the_model_structure%modelboundaries(time)

! get problem matrix and residual vector
     call sparse_smatrix%startAssembly()
     call the_model_structure%modelfk(q1, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q1)

! computing g and dg for the whole model_structure at q0 and q1.
     call the_model_structure%modelgdg(q1, q1, sparse_smatrix, indicesg, indicesq, lambda)
	
! computing constraint forces corresponding to q
     call the_model_structure%modelfL(lambda,sparse_smatrix)
        
     sparse_smatrix%csr_values(sparse_csr_compr_source_kqq) = sparse_smatrix%csr_values(sparse_csr_compr_source_kqq)*2.0d0
     sparse_csr_values_K_stat = 0.0d0
     sparse_csr_values_K_stat = sparse_smatrix%csr_values(sparse_csr_compr_source_K_stat)
        
! Solving linear system of equations
     rhs_unbalance_force = 0.0d0
     rhs_unbalance_force(indicesq) = the_model_structure%fqext

     error = dss_factor_real(handle, MKL_DSS_POSITIVE_DEFINITE, sparse_csr_values_K_stat)
     error = dss_solve_real_d(handle, MKL_DSS_DEFAULTS, rhs_unbalance_force, 1, deltaxvector)
     if (boolean_determinant) then
        statIn = 'determinant'
        call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
        error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
     end if
     
     delta_q_r      = deltaxvector(indicesq)
     delta_lambda_r = deltaxvector(indicesg_stat)
     
     stiffness_parameter = dot_product(delta_q_r,the_model_structure%fqext)/dot_product(delta_q_r,delta_q_r)
     crit_u = dot_product(delta_q_r,D_q)
      
! Our criterium to determine change in (predictor-) load direction
     if (crit_u >= 0.0d0) then
       sign = 1.0d0
     elseif (crit_u < 0.0d0) then
       sign = -1.0d0
     end if

! iterative adjustment of arc-length
     if (settings%number_desired_iteration == 0) then
      DL_arc_n = DL_arc0
     else 
       DL_arc_n = I_n/real(iteration)*DL_arc0
     end if
     
     !print '(3X,A17)', 'predictor-step...'
     !print '(5X,A22,1X,E10.2E2)', 'matrix determinant   =', statOut(1)
     !print '(5X,A22,1X,E10.2E2)', 'stiffnessparameter   =', stiffness_parameter
     !print '(5X,A22,1X,E10.2E2)', 'predictor direction  =', crit_u
     !print '(5X,A22,1X,E10.2E2)', 'predictor arc-length =', DL_arc_n
    
     D_alfa = sign*DL_arc_n/sqrt(dot_product(delta_q_r,delta_q_r))
     D_q    = delta_q_r*D_alfa
     
! additive actualization of generalized displacements
     q1     = qn      + D_q
     lambda = lambdan + delta_lambda_r*D_alfa
     alfa   = alfan   + D_alfa
     
! start Newton iteration
     iteration = 0
     residuum  = 1.0d0
     print '(3X,A17)', 'corrector-step...'
     do while((residuum > settings%tolerance).and.(iteration <= settings%iterationlimit))
        iteration = iteration+1
		
! get problem matrix and residual vector
        call sparse_smatrix%startAssembly()
        call the_model_structure%modelfk(q1, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q1)
		
! computing g and dg for the whole model_structure at q0 and q1.
        call the_model_structure%modelgdg(q1, q1, sparse_smatrix, indicesg, indicesq, lambda)
	
! computing constraint forces corresponding to q
        call the_model_structure%modelfL(lambda,sparse_smatrix)
        
! construction of the residual vector      
        rhs_unbalance_force(indicesq)      = the_model_structure%fqint + the_model_structure%fqL - the_model_structure%fqext*alfa
        rhs_unbalance_force(indicesg_stat) = the_model_structure%g

        sparse_smatrix%csr_values(sparse_csr_compr_source_kqq) = sparse_smatrix%csr_values(sparse_csr_compr_source_kqq)*2.0d0
        sparse_csr_values_K_stat = 0.0d0
        sparse_csr_values_K_stat = sparse_smatrix%csr_values(sparse_csr_compr_source_K_stat)

! Solving linear system of equations
        rhs_solver(:,:)        = 0.0d0
        rhs_solver(:,1)        = -rhs_unbalance_force
        rhs_solver(indicesq,2) = the_model_structure%fqext
        
        error = dss_factor_real(handle, MKL_DSS_POSITIVE_DEFINITE, sparse_csr_values_K_stat)
        error = dss_solve_real_d(handle, MKL_DSS_DEFAULTS, rhs_solver, 2, delta_x_solver)
      
        delta_q_g    = delta_x_solver(indicesq,1)
        delta_q_r    = delta_x_solver(indicesq,2)
        
        delta_lambda_g = delta_x_solver(indicesg_stat,1)
        delta_lambda_r = delta_x_solver(indicesg_stat,2)
        
! arc-length method along a linearized sphere
        if (settings%arc_len_method_flag == 1) then
          ffun       = dot_product(D_q,D_q) + D_alfa*D_alfa - DL_arc_n*DL_arc_n
          f_q        = 2.0d0*D_q
          f_alfa     = 2.0d0*D_alfa
          delta_alfa = - ( ffun + dot_product(f_q,delta_q_g) ) / ( dot_product(f_q,delta_q_r) + f_alfa)

! arc-length method along tangent plane
        else if (settings%arc_len_method_flag == 2) then
          ffun       = 0.0d0
          f_q        = D_q
          f_alfa     = D_alfa
          delta_alfa = - ( ffun + dot_product(f_q,delta_q_g) ) / ( dot_product(f_q,delta_q_r) + f_alfa)
          
! exact arc-length solution
        else if (settings%arc_len_method_flag == 3) then
          a1 = dot_product(delta_q_r,delta_q_r)
          a2 = 2.0d0*dot_product(delta_q_r,D_q + delta_q_g)
          a3 = dot_product(delta_q_g , delta_q_g) + 2.0d0*dot_product(delta_q_g , D_q)
          
          if (a2*a2 < a1*a3) then
            delta_alfa = - a3/a2
          else
            delta_alfa_1   = - ( (a2) + sqrt(a2*a2 - a1*a3) ) / a1
            delta_alfa_2   = - ( (a2) - sqrt(a2*a2 - a1*a3) ) / a1
            a4 = dot_product(D_q, delta_q_g) + dot_product(D_q, D_q)
            a5 = dot_product(D_q, delta_q_r)
            cos1 = (a4 + a5 * delta_alfa_1) / (DL_arc_n*DL_arc_n)
            cos2 = (a4 + a5 * delta_alfa_2) / (DL_arc_n*DL_arc_n)
            delta_alfa = delta_alfa_1
            if     (cos1 >= 0.0d0 .AND. cos2 < 0.0d0) then
              delta_alfa = delta_alfa_1
            elseif (cos1 < 0.0d0 .AND. cos2 >= 0.0d0) then
               delta_alfa = delta_alfa_2
            elseif (cos1 >= 0.0d0 .AND. cos2 >= 0.0d0) then
              if (cos2 > cos1) then
                delta_alfa = delta_alfa_2
              end if
            end if
          end if

        end if
        
        delta_q      = delta_q_g      + delta_alfa*delta_q_r
        delta_lambda = delta_lambda_g + delta_alfa*delta_lambda_r
        
        D_q        = D_q + delta_q
        D_alfa     = D_alfa + delta_alfa          
        
        deltaxvector(indicesq)      = delta_q
        deltaxvector(indicesg_stat) = delta_lambda
        
        residuum = dsqrt(dot_product(deltaxvector, deltaxvector))/dsqrt(dot_product(q1, q1) + dot_product(lambda, lambda))
        if (boolean_output_write .eqv. .TRUE.) then
          print '(5X,A9,1X,E10.2E2,1X,A18,1X,E10.2E2)', 'residuum:', residuum, 'unbalanced forces:',dot_product(rhs_unbalance_force, deltaxvector)
        end if

! additive actualization of generalized displacements and Lagrange multiplier
        q1     = qn      + D_q
        alfa   = alfan   + D_alfa
        lambda = lambda  + delta_lambda
        
        if (isnan(residuum)) then
          the_model_structure%boolean_abort = .True.
          exit
        endif

     end do

! Abort criteria
    if (the_model_structure%boolean_abort) exit
    if (residuum > settings%tolerance .or. isnan(residuum)) then
        if (settings%iterationlimit /= 1) then
          print*, 'solution failed to converge at step ', inos, 'time', time
          the_model_structure%boolean_abort = .True.
        end if 
        exit
    end if

! Actualization of the old soltion
    qn      = q1
    alfan   = alfa
    D_qn    = D_q
    lambdan = lambda
     
! writing in result file
    if (boolean_output_write) then
      call the_model_structure%output_write(settings%step, alfan, alfan - the_model_structure%time_t, inos, settings%i_simutype, iteration, qn, v1, lambdan, the_model_structure%stress)
    end if
    
    print '(3X,A5,1X,F10.3,1X,A11,1X,I2,1X,A9,1X,E10.2E2)','time:',alfan,'iterations:',iteration,'residuum:',residuum

  end do
  
  if (the_model_structure%boolean_abort) then
    goto 999 
  end if
  
! Updating model_structure variables
  the_model_structure%q_t            = qn
  the_model_structure%qdot_t         = v1
  the_model_structure%lambda_t       = lambdan
  the_model_structure%time_t         = time
  
  the_model_structure%kenergy_t      = 0.0d0
  the_model_structure%strainenergy_t = 0.0d0
  the_model_structure%gravenergy_t   = 0.0d0
  the_model_structure%penergy_t      = 0.0d0
  the_model_structure%tenergy_t      = 0.0d0
  
  the_model_structure%lmomentum_t    = 0.0d0
  the_model_structure%amomentum_t    = 0.0d0
  
  the_model_structure%fqext_t      = the_model_structure%fqext
  the_model_structure%fqint_t      = the_model_structure%fqint
  the_model_structure%fqdyn_t      = 0.0d0
  the_model_structure%fqL_t        = the_model_structure%fqL
  the_model_structure%fv_t         = the_model_structure%fqext
  
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
  
! Deallocating  
  deallocate(deltaxvector)
  deallocate(rhs_unbalance_force)
  deallocate(sparse_csr_rowIndices_K_stat)
  deallocate(sparse_csr_columns_K_stat)
  deallocate(sparse_csr_values_K_stat)
  deallocate(sparse_csr_compr_source_K_stat)
  
  return
  end subroutine static_solver_arc