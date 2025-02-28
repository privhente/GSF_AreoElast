subroutine fsi_static_solver(the_model_structure, the_model_aero, the_model_fsi, settings, boolean_output_write)
! Importing classes
  use my_math_structure,     only: cross, eye, outer
  use class_model_structure, only: model_structure, step
  use class_model_aero,      only: model_aero
  use class_model_fsi,       only: model_fsi
  use solver_variables
  use mkl_dss
  use class_sparse_matrix
  use solver_functions
  use my_FileIO_functions
  
  implicit none
  
  type(model_structure), intent(inout) :: the_model_structure
  type(model_aero),      intent(inout) :: the_model_aero
  type(model_fsi),       intent(inout) :: the_model_fsi
  type(step),            intent(in)    :: settings
  type(MKL_DSS_HANDLE)	               :: handle

  real(kind = 8), allocatable :: fqgrav0(:), fqext0(:), statOut(:)
  real(kind = 8), allocatable :: amatrix(:,:)
  
  logical                     :: boolean_sparse = .TRUE., boolean_dfae = .FALSE.
  logical                     :: boolean_determinant = .TRUE.  
  logical, intent(in)         :: boolean_output_write
  logical                     :: boolean_update_aero = .FALSE.
  
  real(kind = 8), dimension(:), allocatable :: fae(:), fae_n(:), fae_n_proj(:), unit_fae(:)
  real(kind = 8)                            :: vec_norm_w
  real(kind = 8)                            :: LF_aero, LF_Struct
  real(kind = 8), allocatable               :: ml_fsi_t(:)
  
  real(kind = 8), dimension(:), allocatable :: csr_values_stat, csr_values_kqq 
  integer, dimension(:), allocatable :: csr_columns_stat, csr_rowIndices_stat, indices_qg, csr_compr_source_stat, csr_columns_kqq , csr_rowIndices_kqq , csr_compr_source_kqq
  real(kind = 8), allocatable :: sparse_csr_values_c(:), sparse_csr_values_kae(:), sparse_smatrix_csr_values(:)
  integer, allocatable        :: sparse_csr_columns_c(:), sparse_csr_columns_kae(:), sparse_smatrix_csr_columns(:)
  integer, allocatable        :: sparse_csr_rowIndices_c(:), sparse_csr_rowIndices_kae(:), sparse_smatrix_csr_rowIndices(:), temp_csr_rowIndices(:)
  integer :: n_rowIndices, nNonZeros
  integer :: convert_job_den_to_csr(8)
  real(kind = 8) :: det_pow, det_base, determinant
  real(kind = 8), allocatable, dimension(:,:) :: dense_matrix_S, dense_matrix_Kae, dense_matrix_Ks
  integer :: convert_job_csr_dense(8)
  
  integer :: UnIn_det, io_error
  logical :: boolean_itsopen
  
  integer(kind = 8) :: sim_clock_rate, delta_sim_clock_start_aero, delta_sim_clock_end_aero, delta_sim_clock_start_structure, delta_sim_clock_end_structure
  real(kind = 8) :: delta_sim_time_aero = 0.0d0, delta_sim_time_structure = 0.0d0
  
!< Dimensioning result vectors for dynamic calculation
  nMatSize     = tncoordinates + trconstraints
  n_rowIndices = the_model_structure%tncoordinates+1

!< allocating arrays for sparse operations
  allocate(sparse_csr_values_c(the_model_aero%tncoordinates*the_model_structure%tncoordinates))
  allocate(sparse_csr_columns_c(the_model_aero%tncoordinates*the_model_structure%tncoordinates))
  allocate(sparse_csr_rowIndices_c(n_rowIndices))
  
  sparse_csr_values_c = 0.0d0
  sparse_csr_columns_c = 0
  sparse_csr_rowIndices_c = 0
  
  allocate(sparse_csr_values_kae(the_model_structure%tncoordinates*the_model_structure%tncoordinates))
  allocate(sparse_csr_columns_kae(the_model_structure%tncoordinates*the_model_structure%tncoordinates))
  allocate(sparse_csr_rowIndices_kae(n_rowIndices))

  sparse_csr_values_kae = 0.0d0
  sparse_csr_columns_kae = 0
  sparse_csr_rowIndices_kae = 0

  allocate(sparse_smatrix_csr_values(nMatSize*nMatSize))
  allocate(sparse_smatrix_csr_columns(nMatSize*nMatSize))
  allocate(sparse_smatrix_csr_rowIndices(nMatSize+1))
  
  sparse_smatrix_csr_values = 0.0d0
  sparse_smatrix_csr_columns = 0
  sparse_smatrix_csr_rowIndices = 0

  allocate(temp_csr_rowIndices(nMatSize+1))
  temp_csr_rowIndices = 0

  allocate(statOut(5))
  statOut = 0.0d0

!< Switching warnings off. PARDISO compatibility problem with new version OPENMP! Maybe resolved in newer version
  !call kmp_set_warnings_off() 
  
  boolean_dfae = .FALSE.
  if (the_model_fsi%flag_dfae .eq. 1) boolean_dfae = .TRUE.

  convert_job_den_to_csr(1) = 0
  convert_job_den_to_csr(2) = 1
  convert_job_den_to_csr(3) = 1
  convert_job_den_to_csr(4) = 2
  convert_job_den_to_csr(5) = the_model_aero%tncoordinates*the_model_aero%tncoordinates  
  convert_job_den_to_csr(6) = 1

  convert_job_csr_dense = 0
  convert_job_csr_dense(1) = 1
  convert_job_csr_dense(2) = 1
  convert_job_csr_dense(3) = 1
  convert_job_csr_dense(4) = 2
  
!< allocating solution vector  
  allocate(deltaxvector(nMatSize))
  allocate(rhs_unbalance_force(nMatSize))
  allocate(fqgrav0(tncoordinates))
  allocate(fqext0(tncoordinates))
  allocate(indices_qg(nMatSize))
  allocate(amatrix(the_model_aero%tnrings,the_model_aero%tnrings))
  allocate(fae(tncoordinates))
  allocate(fae_n(tncoordinates))
  allocate(fae_n_proj(tncoordinates))
  allocate(unit_fae(tncoordinates))
  allocate(ml_fsi_t(the_model_structure%nnodes12*6+the_model_structure%nnodes6*3))
  
  indices_qg = 0.0d0
  indices_qg = reshape((/indicesq,indicesg/),(/nMatSize/))

!< Creating and extracting sparse matrix from dynamic system matrix for static analysis
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, csr_values_stat, csr_columns_stat, csr_rowIndices_stat, csr_compr_source_stat, indices_qg, indices_qg)
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, csr_values_kqq , csr_columns_kqq , csr_rowIndices_kqq , csr_compr_source_kqq , indicesq,   indicesq)  
  
!< Initialization of field variables
  q0                  = the_model_structure%q_t
  lambda0             = the_model_structure%lambda_t
  fqext0              = the_model_structure%fqext0_t
  fqgrav0             = 0.0d0
  fae_n               = 0.0d0
  fae                 = 0.0d0
  fae_n_proj          = 0.0d0
  unit_fae            = 0.0d0
  deltaxvector        = 0.0d0
  rhs_unbalance_force = 0.0d0
  amatrix             = 0.0d0
  ml_fsi_t            = 0.0d0
  
  the_model_structure%fqgrav(:) = 0.0d0
  if (settings%gravityflag == 1) then
    call the_model_structure%modelfgrav()
  end if
  fqgrav0 = the_model_structure%fqgrav

!< settings for time stepping in structural and aero time steps  
  the_model_structure%deltat = settings%deltat
  if (the_model_structure%deltat .ge. the_model_aero%deltat) the_model_structure%deltat = the_model_aero%deltat

!< set airflow flag to TRUE
  the_model_aero%airflow%boolean_actualize = .TRUE.
  
!< set unsteady flag to false to switch off the unsteady term in force calculation
  the_model_aero%boolean_unsteady_term = .FALSE.

!< if not strong fsi then also not linearization of aerodynamic loads
  if (.not. the_model_fsi%boolean_strongfsi) boolean_dfae = .FALSE.
  
!< for file in case determinant shall be written  
  call GetNewUnit(UnIn_det)
  
!< ------------- START AERODYNAMIC TIME STEP --------------------------------------------------------------------------------------------------------
  do while (the_model_aero%tcounter+1 .le. the_model_aero%nsteps)
    
    call system_clock(delta_sim_clock_start_aero, sim_clock_rate)
    
    !< update wake and tcounter
    if (the_model_aero%tcounter .gt. 0) then
      call the_model_aero%update_wake()
    end if
    
    the_model_aero%time = the_model_aero%tcounter*the_model_aero%deltat
    LF_aero = the_model_fsi%calc_load_factor(the_model_aero%time)
    
    !< set initial free-field velocity for external flow field; for prescirbed flow-field from file, this is valid for outside of flow-field grid
    the_model_aero%airflow%vinfinity = the_model_aero%airflow%dir * the_model_aero%evaluateamplitude(the_model_aero%time)
    
    !< Starting loop over structural time steps
    the_model_structure%time_t             = the_model_aero%time - the_model_aero%deltat
    the_model_structure%settings(1)%totalt = the_model_aero%time
    time                                   = the_model_aero%time - the_model_aero%deltat 
    
!< ------------- START STRUCTURAL TIME STEP --------------------------------------------------------------------------------------------------------
    inos  = IDNINT(the_model_aero%deltat/the_model_structure%deltat)
    istep = 0
    delta_sim_time_structure = 0.0d0
    print '(2X,A10,1X,I4,1X,A5,1X,F10.3,1X,A11,I3)', 'Aero step:', the_model_aero%tcounter, 'time:', the_model_aero%time, 'load steps:', inos
    do while (istep .lt. inos)
      
      !< time measure of each structural time step
      call system_clock(delta_sim_clock_start_structure, sim_clock_rate)
      
      istep   = istep + 1
      time    = the_model_structure%time_t + istep*the_model_structure%deltat
      
      ! load factor for structural load stepping
      LF_Struct = dble(istep)/dble(inos)
            
      !< set boolean_update_aero to true for to calculate aerodyn in first step
      boolean_update_aero = .TRUE.

    !< Get model boundaries at current time increment
      if (time .le. 0.0d0) then
        call the_model_structure%modelloads(0.0d0, ml_t, sl_t)
        call the_model_structure%modelboundaries(0.0d0)
        time = 0.0d0
        istep = inos
      else
        call the_model_structure%modelloads(time, ml_t, sl_t)
        call the_model_structure%modelboundaries(time)
      end if
      
      !< gravity forces
      if (settings%gravityflag .eq. 1) then
        the_model_structure%fqgrav = fqgrav0
      end if
      
      !< Configuration for positions and velocities at current time step 
      q1     = q0
      lambda = lambda0
      
!< ------------- START NEWTON ITERATION --------------------------------------------------------------------------------------------------------
      iteration = 0
      residuum  = 1.0d0
      do while((residuum > settings%tolerance).and.(iteration < settings%iterationlimit))
        
        ! check simulation and iteration time
        call system_clock(delta_sim_clock_start_structure, sim_clock_rate)
        
        iteration = iteration+1

        !< computing structural vectors and matrices
        call sparse_smatrix%startAssembly()
        call the_model_structure%modelfk(q1, ml_fsi_t, sl_t, sparse_smatrix, indicesq, indicesv, q1)
        call the_model_structure%modelgdg(q1, q1, sparse_smatrix, indicesg, indicesq, lambda, q0)
        call the_model_structure%modelfL(lambda,sparse_smatrix)
        
        !< multiply the sparse entries by 2.0d0, because the static solution is at tn+1, not tn05
        sparse_smatrix%csr_values(csr_compr_source_kqq) = sparse_smatrix%csr_values(csr_compr_source_kqq)*2.0d0
        csr_values_stat = 0.0d0
        csr_values_stat = sparse_smatrix%csr_values(csr_compr_source_stat)
        
        if (boolean_update_aero) then
          !< computing aerodynamical forces and tangent matrix (if switched on)
          
          !< mapping of vortex-lattice kinematics according to deformed FE-mesh
          call the_model_fsi%actualize(the_model_aero, q1, the_model_structure%tncoordinates) ! assign structure's geometry to aerogrid's geometry
          
          !< update vortex-lattice (vortex rings and segments) kinematics with mapping coming from the_model_fsi%actualize
          call the_model_aero%actualize(.TRUE., .FALSE.) ! actualize geometry
          
          ! update connecting wake first row segments and first row rings clued to separation edge
          call the_model_aero%wake_update_first_row(.True.,.False.)
         
          !< update the external flow field only for the first structural step and first Newton iteration
          if (istep .eq. 1 .and. iteration .eq. 2) the_model_aero%airflow%boolean_actualize = .FALSE.
        
          !< evalute aerodynamics (uvlm)
          call the_model_aero%evaluaterhs(the_model_aero%cutoff)
          call the_model_aero%evaluateamatrix(the_model_aero%cutoff)
          the_model_aero%circulations = the_model_aero%rhs
          amatrix = the_model_aero%amatrix
          call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
          call the_model_aero%actualize(.FALSE., .TRUE.)  ! actualize ring circulations
        
          !< Determine aerodynamical forces and tangent matrix 
          call the_model_aero%modelfae()
        
          !< transform aerodynamic forces at aerodynamic nodes to structural nodes
          call mkl_dcsrmv('T', the_model_aero%tncoordinates, the_model_structure%tncoordinates, 1.0d0, 'G  F  ', &
            the_model_fsi%sparse_Tn%csr_values, the_model_fsi%sparse_Tn%csr_columns, &
            the_model_fsi%sparse_Tn%csr_rowIndices(1:size(the_model_fsi%sparse_Tn%csr_rowIndices) - 1), &
            the_model_fsi%sparse_Tn%csr_rowIndices(2:size(the_model_fsi%sparse_Tn%csr_rowIndices)), &
            the_model_aero%faen, 0.0d0, the_model_fsi%fae)

!< ToDo>: HERE CHANGES REGARDING LOADS WHILE WEAK COUPLING:
          !< initialize the vector of aerodynamic material amplitude loads
          if (.not. the_model_fsi%boolean_strongfsi) then
            !< calculate vector of amplitudes as input for material loads: the_model_fsi%ml_aero_t
            call the_model_fsi%amplitude_material_loads(the_model_structure, q1)
            ml_fsi_t = ml_t + the_model_fsi%ml_aero_t
          end if
!< UNTIL HERE!
      
          !< perform lineatization of aerodynamic forces/moments
          !< extend structural sparse matrix by sparse aerodynamic tangent matrix
          if (boolean_dfae) then
            call the_model_aero%modeldfae(boolean_dfae)

            !< write dfae in sparse format
            if (allocated(the_model_aero%sparse_dfae%csr_values)) deallocate(the_model_aero%sparse_dfae%csr_values)
            if (allocated(the_model_aero%sparse_dfae%csr_columns)) deallocate(the_model_aero%sparse_dfae%csr_columns)
            if (allocated(the_model_aero%sparse_dfae%csr_rowIndices)) deallocate(the_model_aero%sparse_dfae%csr_rowIndices)
          
            nNonZeros = count(the_model_aero%dfae(1:the_model_aero%tncoordinates,1:the_model_aero%tncoordinates) .ne. 0.0d0)
            allocate(the_model_aero%sparse_dfae%csr_values(nNonZeros))
            allocate(the_model_aero%sparse_dfae%csr_columns(nNonZeros))
            allocate(the_model_aero%sparse_dfae%csr_rowIndices(the_model_aero%tncoordinates+1))
        
            call mkl_ddnscsr(convert_job_den_to_csr, the_model_aero%tncoordinates, the_model_aero%tncoordinates, & 
                                  the_model_aero%dfae(1:the_model_aero%tncoordinates,1:the_model_aero%tncoordinates), the_model_aero%tncoordinates, & 
                                  the_model_aero%sparse_dfae%csr_values, the_model_aero%sparse_dfae%csr_columns, the_model_aero%sparse_dfae%csr_rowIndices, info)
        
            !< Performing sparse multiplication to get aerodynamical tangent matrix with repect to structural coordinates
            !< [C] = Tn' * dfae_n; (q x n)
            sparse_csr_values_c  = 0.0d0
            sparse_csr_columns_c = 0
            call mkl_dcsrmultcsr('T', 0, 8, the_model_aero%tncoordinates, the_model_structure%tncoordinates, the_model_aero%tncoordinates, &
                                            the_model_fsi%sparse_Tn%csr_values, the_model_fsi%sparse_Tn%csr_columns, the_model_fsi%sparse_Tn%csr_rowIndices, &
                                            the_model_aero%sparse_dfae%csr_values, the_model_aero%sparse_dfae%csr_columns, the_model_aero%sparse_dfae%csr_rowIndices, &
                                            sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, the_model_aero%tncoordinates*the_model_structure%tncoordinates, info)
        
            !< [kqq,kqv] = C*T2n; (q x q)
            sparse_csr_values_kae  = 0.0d0
            sparse_csr_columns_kae = 0
            nNonZeros              = sparse_csr_rowIndices_c(n_rowIndices)-1
            call mkl_dcsrmultcsr('N', 0, 8, the_model_structure%tncoordinates, the_model_aero%tncoordinates, the_model_structure%tncoordinates, &
                                            sparse_csr_values_c(1:nNonZeros), sparse_csr_columns_c(1:nNonZeros), sparse_csr_rowIndices_c, &
                                            the_model_fsi%sparse_Tn%csr_values, the_model_fsi%sparse_Tn%csr_columns, the_model_fsi%sparse_Tn%csr_rowIndices, &
                                            sparse_csr_values_kae, sparse_csr_columns_kae, sparse_csr_rowIndices_kae, the_model_structure%tncoordinates*the_model_structure%tncoordinates, info)
        
            !< Merging aerodynamical and structural tangent matrices
            temp_csr_rowIndices                 = sparse_csr_rowIndices_kae(n_rowIndices)       !< adjust the dimension of aerodynamic matrix to dimension of system matrix, i.e. matSize
            temp_csr_rowIndices(1:n_rowIndices) = sparse_csr_rowIndices_kae
            sparse_smatrix_csr_values           = 0.0d0
            sparse_smatrix_csr_columns          = 0
            nNonZeros                           = sparse_csr_rowIndices_kae(n_rowIndices)-1
            call mkl_dcsradd('N',0,4,nMatSize,nMatSize, &
                                           csr_values_stat, csr_columns_stat, csr_rowIndices_stat, & 
                                          -1.0d0*LF_Aero*LF_Struct,sparse_csr_values_kae(1:nNonZeros),sparse_csr_columns_kae(1:nNonZeros),temp_csr_rowIndices, &
                                           sparse_smatrix_csr_values,sparse_smatrix_csr_columns,sparse_smatrix_csr_rowIndices,nMatSize*nMatSize,info)
          else
            sparse_smatrix_csr_rowIndices = csr_rowIndices_stat
            sparse_smatrix_csr_columns(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1) = csr_columns_stat
            sparse_smatrix_csr_values(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1)  = csr_values_stat
          end if
          
          if (inos .gt. 1) then 
              the_model_fsi%fae = the_model_structure%fqext
              the_model_structure%fqext = 0.0d0
          end if
          
          !< vector of generalized aerodynamic forces for current time step - TO DO - UPDATE THE LINEARIZATION FOR THIS LOAD STEPPING
          if (.not. boolean_dfae) then
            fae = fae_n_proj + (the_model_fsi%fae - fae_n_proj)*LF_Struct
          else !< just a simplification here, it has to be addjusted
            fae = fae_n + (the_model_fsi%fae - fae_n)*LF_Struct
          end if
            
          !< in case of weak coupling, switch off boolean_strongfsi
          if (.not. the_model_fsi%boolean_strongfsi) boolean_update_aero = .FALSE.

        else
          sparse_smatrix_csr_rowIndices = csr_rowIndices_stat
          sparse_smatrix_csr_columns(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1) = csr_columns_stat
          sparse_smatrix_csr_values(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1)  = csr_values_stat
        end if
        
        if (.not. the_model_fsi%boolean_strongfsi) then
            the_model_fsi%fae = the_model_structure%fqext
            the_model_structure%fqext = 0.0d0
        end if
          
        !< construction of the residual vector
        rhs_unbalance_force(indicesq)      =  the_model_structure%fqint + the_model_structure%fqdyn + the_model_structure%fqL - & 
                                              the_model_structure%fqgrav - the_model_structure%fqext - & 
                                              LF_Aero*fae
        
        rhs_unbalance_force(indicesg_stat) = the_model_structure%g
        
        !< Solving linear system of equation
        nNonZeros = sparse_smatrix_csr_rowIndices(nMatSize+1)-1
        if (boolean_sparse) then
          error = dss_create(handle, MKL_DSS_DEFAULTS)
          error = dss_define_structure(handle, MKL_DSS_NON_SYMMETRIC, sparse_smatrix_csr_rowIndices, nMatSize, nMatSize, sparse_smatrix_csr_columns(1:nNonZeros), nNonZeros )
          error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
        end if
        
        call solver(nMatSize, handle, sparse_smatrix_csr_values(1:nNonZeros), sparse_smatrix_csr_columns(1:nNonZeros), sparse_smatrix_csr_rowIndices, rhs_unbalance_force, 1, deltaxvector, boolean_sparse, statOut, boolean_determinant)
        
        if (boolean_sparse) then
          error = dss_delete(handle, MKL_DSS_MSG_LVL_WARNING+MKL_DSS_TERM_LVL_ERROR)
        end if
        
        !< print determinat
        if (boolean_determinant) then
          inquire(unit=UnIn_det, opened=boolean_itsopen)
          if (.not. boolean_itsopen) then
            call OpenFOutFileToRead(UnIn_det, 'determinant.dres', io_error)
          end if
          write(UnIn_det,'(1000000f30.15)') statOut(1), statOut(2)
        end if
    
        if (the_model_aero%time .eq. 0.0d0) then
          deltaxvector = 0.0d0
          residuum     = 0.0d0
          exit
        end if
        
        residuum = dsqrt(dot_product(deltaxvector, deltaxvector))/dsqrt(dot_product(q1, q1)+dot_product(lambda, lambda))
        print '(8X,A9,1X,E10.2E2,1X,A17,1X,E10.2E3)', 'residuum:', residuum, 'convergence work:',abs(dot_product(rhs_unbalance_force, deltaxvector))
        
        if (isnan(residuum)) then
          residuum = 1.0d0
          the_model_structure%boolean_abort = .True.
          exit
        endif
        
        !< additive actualization of generalized displacements and Lagrange multiplier
        q1     = q1     + deltaxvector(indicesq)
        lambda = lambda + deltaxvector(indicesg_stat)

        ! check simulation and iteration time
        call system_clock(delta_sim_clock_end_structure)
        delta_sim_time_structure = delta_sim_time_structure + (real(delta_sim_clock_end_structure -delta_sim_clock_start_structure))/real(sim_clock_rate)
      end do
     
      !< abort criteria
      if (residuum > settings%tolerance) then
        print*, '... solution failed to converge at time', time
        if (boolean_dfae) then
          the_model_structure%boolean_abort = .TRUE.
          exit
        else
          !print*, '... restart load step with linearized aerodynamic forces'
          boolean_dfae = .TRUE.
          istep = istep - 1
          the_model_structure%boolean_abort = .True.
          exit
        end if
      else
        !< actualization of the old soltion     
        q0      = q1
        lambda0 = lambda
        print '(5X,A22,1X,F10.5,1X,A11,1X,I2,1X,A9,1X,E10.2E3,1X,A3,F10.5)','converged solution at:',time,'iterations:',iteration,'residuum:',residuum,'t:',delta_sim_time_structure        
      end if
      
    end do
    
    if (the_model_structure%boolean_abort) then
      exit
    end if
    
    !< update aero grid to current solution
    call the_model_fsi%actualize(the_model_aero, q1, the_model_structure%tncoordinates) ! assign structure's geometry to aerogrid's geometry
    call the_model_aero%actualize(.TRUE., .FALSE.) ! actualize geometry
    call the_model_aero%wake_update_first_row()
    
    !< update wake data for current aero time step
    the_model_aero%airflow%boolean_actualize = .TRUE. 
    call the_model_aero%wake_velocity(the_model_aero%cutoff)
    the_model_aero%oldcirculations = the_model_aero%circulations

    !< calculate invariants of time current step
    call the_model_structure%modelinvariants(q1, v1)

    !< write results to file
    call the_model_aero%output_write()
    call the_model_structure%output_write(the_model_aero%tcounter,the_model_aero%time, the_model_aero%time, the_model_aero%tcounter, settings%i_simutype, iteration, q1, v1, lambda, the_model_structure%stress_t)
    
    !< update wake and tcounter
    the_model_aero%tcounter = the_model_aero%tcounter+1
    
    !< update fae_n
    fae_n = fae
    
    call system_clock(delta_sim_clock_end_aero)
    delta_sim_time_aero = (real(delta_sim_clock_end_aero - delta_sim_clock_start_aero))/real(sim_clock_rate)
    print '(2X,A31,1X,F10.5,1X,A3,F10.5)','... Aero converged solution at:',the_model_aero%time,'t:',delta_sim_time_aero
    
  end do
  
  if (boolean_determinant) then
    close(UnIn_det)
  end if
  
  if (the_model_structure%boolean_abort) then
    print '(5X,A23,1X,F10.5,1X,A11,1X,I2,1X,A9,1X,E10.2E3)','no convergence at time:',time,'iterations:',iteration,'residuum:',residuum
    go to 999
  end if
  
  !< Updating model variables
  the_model_structure%q_t            = q1 
  the_model_structure%lambda_t       = lambda
  the_model_structure%time_t         = time
  
  the_model_structure%kenergy_t      = the_model_structure%kenergy
  the_model_structure%strainenergy_t = the_model_structure%strainenergy
  the_model_structure%gravenergy_t   = the_model_structure%gravenergy
  the_model_structure%penergy_t      = the_model_structure%penergy
  the_model_structure%tenergy_t      = the_model_structure%tenergy
  
  the_model_structure%lmomentum_t    = the_model_structure%lmomentum
  the_model_structure%amomentum_t    = the_model_structure%amomentum
  
  the_model_structure%fqext_t      = the_model_structure%fqext
  the_model_structure%fqext0_t     = fqext0
  the_model_structure%fqint_t      = the_model_structure%fqint
  the_model_structure%fqdyn_t      = the_model_structure%fqdyn
  the_model_structure%fqL_t        = the_model_structure%fqL
  the_model_structure%fv_t         = the_model_structure%fv

  if (settings%output_flag == 1) then
    print *, "Output: Writing matrices..."
    
    call writeIntVectorToFile(indicesq, 9999,'indicesq.dres')
    call writeIntVectorToFile(indicesg_stat, 9999,'indicesg.dres')
      
    call writeIntVectorToFile(csr_columns_stat(1:csr_rowIndices_stat(nMatSize+1)-1),   9999,'csr_columns_Kstat.dres')
    call writeIntVectorToFile(csr_rowIndices_stat,9999,'csr_rowIndices_Kstat.dres')
    call writeRealVectorToFile(csr_values_stat(1:csr_rowIndices_stat(nMatSize+1)-1),   9999,'csr_values_Kstat.dres')
    
	  call writeIntVectorToFile(sparse_csr_columns_kae(1:sparse_csr_rowIndices_kae(n_rowIndices)-1),   9999,'csr_columns_Kae.dres')
    call writeIntVectorToFile(sparse_csr_rowIndices_kae,9999,'csr_rowIndices_Kae.dres')
    call writeRealVectorToFile(sparse_csr_values_kae(1:sparse_csr_rowIndices_kae(n_rowIndices)-1),9999,'csr_values_Kae.dres')
    
	  call writeIntVectorToFile(sparse_smatrix_csr_columns(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1),   9999,'csr_columns_S.dres')
    call writeIntVectorToFile(sparse_smatrix_csr_rowIndices,9999,'csr_rowIndices_S.dres')
    call writeRealVectorToFile(sparse_smatrix_csr_values(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1),9999,'csr_values_S.dres')
  end if  

999 continue
    
  deallocate(amatrix)
  deallocate(deltaxvector)
  deallocate(indices_qg)
  deallocate(rhs_unbalance_force)
  deallocate(csr_rowIndices_stat)
  deallocate(csr_columns_stat)
  deallocate(csr_values_stat)
  deallocate(csr_compr_source_stat)
  deallocate(fqgrav0)
  deallocate(fqext0)
  return
end subroutine fsi_static_solver
