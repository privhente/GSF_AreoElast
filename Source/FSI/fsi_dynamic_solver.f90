#if Controller_ADDON
subroutine fsi_dynamic_solver(the_model_structure, the_model_aero, the_model_controller, the_model_fsi, settings, boolean_output_write)
#else
subroutine fsi_dynamic_solver(the_model_structure, the_model_aero, the_model_fsi, settings, boolean_output_write)
#endif
  use my_math_structure,      only: cross, eye, outer
  use class_model_structure,  only: model_structure, step
  use class_model_aero,       only: model_aero
  use class_model_fsi,        only: model_fsi
  #if Controller_ADDON
    use class_controller,       only: model_controller
  #endif
  use solver_variables
  use mkl_dss
  use class_sparse_matrix
  use solver_functions
  use my_FileIO_functions
  use my_fsi_functions, only: det
 
  implicit none

  #if Controller_ADDON
    type(model_controller)               :: the_model_controller
  #endif
  type(model_structure), intent(inout) :: the_model_structure
  type(model_aero),      intent(inout) :: the_model_aero
  type(model_fsi),       intent(inout) :: the_model_fsi
  type(step),            intent(in)    :: settings
  type(MKL_DSS_HANDLE)	               :: handle

  real(kind = 8), allocatable :: fqext0(:), statOut(:), dense(:,:), fae(:), fae_n(:), fae_n_proj(:), unit_fae(:)
  real(kind = 8)              :: vec_norm_w
  real(kind = 8)              :: LF_aero, LF_Struct
  real(kind = 8), allocatable :: amatrix(:,:), q05(:), v05(:), v0_aero(:), q0_test(:)

  real(kind = 8), allocatable :: sparse_csr_values_c(:), sparse_csr_values_kae(:), sparse_smatrix_csr_values(:)
  integer, allocatable        :: sparse_csr_columns_c(:), sparse_csr_columns_kae(:), sparse_smatrix_csr_columns(:)
  integer, allocatable        :: sparse_csr_rowIndices_c(:), sparse_csr_rowIndices_kae(:), sparse_smatrix_csr_rowIndices(:), temp_csr_rowIndices(:)
  integer :: n_rowIndices, nNonZeros
  integer :: convert_job_den_to_csr(8), convert_job_csr_to_den(8)

  integer(kind = 8) :: sim_clock_rate, delta_sim_clock_start_aero, delta_sim_clock_end_aero, delta_sim_clock_start_structure, delta_sim_clock_end_structure
  real(kind = 8) :: delta_sim_time_aero = 0.0d0, delta_sim_time_structure = 0.0d0

  integer :: UnIn_det, io_error

  real(kind = 8), allocatable, dimension(:,:) :: dense_dG

  logical, intent(in)         :: boolean_output_write
  logical                     :: boolean_itsopen
  logical                     :: boolean_sparse = .TRUE., boolean_dfae = .FALSE., boolean_determinant = .FALSE.
  logical                     :: boolean_update_aero = .FALSE.
  
!< Switching warnings off. PARDISO compatibility problem with new version OPENMP! Maybe resolved in newer version
  !call kmp_set_warnings_off()

  boolean_dfae = .FALSE.
  if (the_model_fsi%flag_dfae .eq. 1) boolean_dfae = .TRUE.

!< set conversion setting for sparse-dense transformation and vice versa
  convert_job_den_to_csr(1) = 0
  convert_job_den_to_csr(2) = 1
  convert_job_den_to_csr(3) = 1
  convert_job_den_to_csr(4) = 2
  convert_job_den_to_csr(5) = 2*the_model_aero%tncoordinates*the_model_aero%tncoordinates
  convert_job_den_to_csr(6) = 1

  convert_job_csr_to_den(:) = 0
  convert_job_csr_to_den(1) = 1
  convert_job_csr_to_den(2) = 1
  convert_job_csr_to_den(3) = 1
  convert_job_csr_to_den(4) = 2

!< Dimensioning result vectors for dynamic calculation
  nMatSize     = 2*the_model_structure%tncoordinates + trconstraints
  n_rowIndices = the_model_structure%tncoordinates+1

!< allocating arrays for sparse operations
  allocate(sparse_csr_values_c(2*the_model_aero%tncoordinates*the_model_structure%tncoordinates))
  allocate(sparse_csr_columns_c(2*the_model_aero%tncoordinates*the_model_structure%tncoordinates))
  allocate(sparse_csr_rowIndices_c(n_rowIndices))

  sparse_csr_values_c = 0.0d0
  sparse_csr_columns_c = 0
  sparse_csr_rowIndices_c = 0
    
  allocate(sparse_csr_values_kae(2*the_model_structure%tncoordinates*the_model_structure%tncoordinates))
  allocate(sparse_csr_columns_kae(2*the_model_structure%tncoordinates*the_model_structure%tncoordinates))
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
  
  allocate(dense_dG(the_model_aero%tnrings,2*the_model_structure%tncoordinates))

!< allocate arrays for solving governing equations
  allocate(deltaxvector(nMatSize))
  allocate(rhs_unbalance_force(nMatSize))
  allocate(fqext0(tncoordinates))
  allocate(amatrix(the_model_aero%tnrings,the_model_aero%tnrings))
  allocate(q05(tncoordinates))
  allocate(q0_test(tncoordinates))
  allocate(v05(tncoordinates))
  allocate(v0_aero(tncoordinates))
  allocate(fae(tncoordinates))
  allocate(fae_n(tncoordinates))
  allocate(fae_n_proj(tncoordinates))
  allocate(unit_fae(tncoordinates))

!< Initialization of field variables
  q0                  = the_model_structure%q_t
  v0                  = the_model_structure%qdot_t
  lambda0             = the_model_structure%lambda_t
  fqext0              = the_model_structure%fqext0_t
  fae_n               = 0.0d0
  fae                 = 0.0d0
  fae_n_proj          = 0.0d0
  unit_fae            = 0.0d0
  deltaxvector        = 0.0d0
  rhs_unbalance_force = 0.0d0
  amatrix             = 0.0d0
  q05                 = 0.0d0
  v05                 = 0.0d0
  v0_aero             = 0.0d0
  q0_test             = 0.0d0
  
!< gravitiy force vector
  the_model_structure%fqgrav(:) = 0.0d0
  if (settings%gravityflag .eq. 1) then
    call the_model_structure%modelfgrav()
  end if

!< settings for time stepping in structural and aero time steps
  the_model_structure%deltat = settings%deltat
  if (the_model_structure%deltat .ge. the_model_aero%deltat) the_model_structure%deltat = the_model_aero%deltat

!< set airflow flag to TRUE
  the_model_aero%airflow%boolean_actualize = .TRUE.

!< set unsteady flag to false to switch off the unsteady term in force calculation
  the_model_aero%boolean_unsteady_term = .TRUE.

!< if not strong fsi then also not linearization of aerodynamic loads
  if (.not. the_model_fsi%boolean_strongfsi) boolean_dfae = .FALSE.

!< for file in case determinant shall be written
  call GetNewUnit(UnIn_det)

!< ------------- START AERODYNAMIC TIME STEP --------------------------------------------------------------------------------------------------------
  do while (the_model_aero%time .le. the_model_aero%totalt)
    call system_clock(delta_sim_clock_start_aero, sim_clock_rate)

    !< update wake and tcounter
    if (the_model_aero%tcounter .gt. 0) then
      call the_model_aero%update_wake()
    end if

    !< update aero time and load factor for Aero
    the_model_aero%time = the_model_aero%time + the_model_aero%deltat
    LF_aero = the_model_fsi%calc_load_factor(the_model_aero%time)

    !< set initial free-field velocity for external flow field; for prescirbed flow-field from file, this is valid for outside of flow-field grid
    the_model_aero%airflow%vinfinity = the_model_aero%airflow%dir * the_model_aero%evaluateamplitude(the_model_aero%time)

    !< updating structural time step
    the_model_structure%time_t             = the_model_aero%time - the_model_aero%deltat
    the_model_structure%settings(1)%totalt = the_model_aero%time
    time                                   = the_model_aero%time - the_model_aero%deltat

    #if Controller_ADDON
        !< assign controller data to model for next time step
        if (the_model_controller%boolean_abort .eqv. .FALSE.) then
        call the_model_controller%assign_servo_to_model(the_model_structure, q1, q0_test, the_model_aero%time, the_model_aero%deltat)
        end if
    #endif

    !< saving old solution to compute accelerations for required in update_controller
    v0_aero = v0

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
      
      !< set boolean_update_aero to true to calculate aerodynamic loads. In case weak is activated, then aerodynamic loads 
      !< are updated at the beginning of each structural load step
      boolean_update_aero = .TRUE.
      
      !< load factor for structural load stepping
      LF_Struct = dble(istep)/dble(inos)

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

      !< Configuration for positions and velocities at current time step
      q1     = q0
      v1     = v0
      lambda = lambda0

!< ------------- START NEWTON ITERATION --------------------------------------------------------------------------------------------------------
      iteration = 0
      residuum  = 1.0d0
      do while((residuum > settings%tolerance).and.(iteration < settings%iterationlimit))

        !< check simulation and iteration time
        call system_clock(delta_sim_clock_start_structure, sim_clock_rate)

        iteration = iteration+1

        !< computing structural vectors and matrices
        call sparse_smatrix%startAssembly()
        call the_model_structure%modelfk(q0, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q1,  v0, v1)
        call the_model_structure%modelgdg(q1, q0, sparse_smatrix, indicesg, indicesq, lambda)
        call the_model_structure%modelfL(lambda,sparse_smatrix)
        
        !< computing aerodynamical forces and tangent matrix (if switched on)
        if (boolean_update_aero) then
          q05 = 0.5d0 * (q1 + q0)
          v05 = 0.5d0 * (v1 + v0)
          
          !< mapping of vortex-lattice kinematics according to deformed FE-mesh
          call the_model_fsi%actualize(the_model_aero, q05, the_model_structure%tncoordinates, v05) ! assign structure's geometry to aerogrid's geometry
          
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

          !< determine aerodynamical forces
          call the_model_aero%modelfae()

          !< transform aerodynamic forces at aerodynamic nodes to structural nodes
          call mkl_dcsrmv('T', the_model_aero%tncoordinates, the_model_structure%tncoordinates, 1.0d0, 'G  F  ', &
            the_model_fsi%sparse_Tn%csr_values, the_model_fsi%sparse_Tn%csr_columns, &
            the_model_fsi%sparse_Tn%csr_rowIndices(1:size(the_model_fsi%sparse_Tn%csr_rowIndices) - 1), &
            the_model_fsi%sparse_Tn%csr_rowIndices(2:size(the_model_fsi%sparse_Tn%csr_rowIndices)), &
            the_model_aero%faen, 0.0d0, the_model_fsi%fae)
          
          !< perform lineatization of aerodynamic forces/moments and extend structural sparse matrix by sparse aerodynamic tangent matrix
          if (boolean_dfae) then
            !< compute tangential matrices of aerodynamic loads with respect to nodes of boundary element mesh (aero-grid)
            call the_model_aero%modeldfae(boolean_dfae)

            if (allocated(the_model_aero%sparse_dfae%csr_values)) deallocate(the_model_aero%sparse_dfae%csr_values)
            if (allocated(the_model_aero%sparse_dfae%csr_columns)) deallocate(the_model_aero%sparse_dfae%csr_columns)
            if (allocated(the_model_aero%sparse_dfae%csr_rowIndices)) deallocate(the_model_aero%sparse_dfae%csr_rowIndices)

            nNonZeros = count(the_model_aero%dfae .ne. 0.0d0)
            allocate(the_model_aero%sparse_dfae%csr_values(nNonZeros))
            allocate(the_model_aero%sparse_dfae%csr_columns(nNonZeros))
            allocate(the_model_aero%sparse_dfae%csr_rowIndices(the_model_aero%tncoordinates+1))
            convert_job_den_to_csr(5) = nNonZeros
            
            !< converting from dense to sparse format, although matrix is highly dense
            call mkl_ddnscsr(convert_job_den_to_csr, the_model_aero%tncoordinates, 2*the_model_aero%tncoordinates, the_model_aero%dfae, the_model_aero%tncoordinates, &
                                  the_model_aero%sparse_dfae%csr_values, the_model_aero%sparse_dfae%csr_columns, the_model_aero%sparse_dfae%csr_rowIndices, info)

            !< performing sparse multiplication to get aerodynamical tangent matrix with respect to structural coordinates
            !< [C] = Tn' * dfae_n; (q x 2n)
            sparse_csr_values_c     = 0.0d0
            sparse_csr_columns_c    = 0
            sparse_csr_rowIndices_c = 0
            call mkl_dcsrmultcsr('T', 0, 8, the_model_aero%tncoordinates, the_model_structure%tncoordinates, 2*the_model_aero%tncoordinates, &
                                            the_model_fsi%sparse_Tn%csr_values, the_model_fsi%sparse_Tn%csr_columns, the_model_fsi%sparse_Tn%csr_rowIndices, &
                                            the_model_aero%sparse_dfae%csr_values, the_model_aero%sparse_dfae%csr_columns, the_model_aero%sparse_dfae%csr_rowIndices, &
                                            sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, 2*the_model_aero%tncoordinates*the_model_structure%tncoordinates, info)

            !< [kqq,kqv] = C*T2n; (q x 2q)
            sparse_csr_values_kae     = 0.0d0
            sparse_csr_columns_kae    = 0
            sparse_csr_rowIndices_kae = 0
            nNonZeros              = sparse_csr_rowIndices_c(n_rowIndices)-1
            call mkl_dcsrmultcsr('N', 0, 8, the_model_structure%tncoordinates, 2*the_model_aero%tncoordinates, 2*the_model_structure%tncoordinates, &
                                            sparse_csr_values_c(1:nNonZeros), sparse_csr_columns_c(1:nNonZeros), sparse_csr_rowIndices_c, &
                                            the_model_fsi%sparse_T2n%csr_values, the_model_fsi%sparse_T2n%csr_columns, the_model_fsi%sparse_T2n%csr_rowIndices, &
                                            sparse_csr_values_kae, sparse_csr_columns_kae, sparse_csr_rowIndices_kae, 2*the_model_structure%tncoordinates*the_model_structure%tncoordinates, info)

            !< merging aerodynamical and structural tangent matrices by adding both
            temp_csr_rowIndices                 = sparse_csr_rowIndices_kae(n_rowIndices)      !< adjust the dimension of aerodynamic matrix to dimension of system matrix, i.e. matSize
            temp_csr_rowIndices(1:n_rowIndices) = sparse_csr_rowIndices_kae
            sparse_smatrix_csr_values           = 0.0d0
            sparse_smatrix_csr_columns          = 0
            nNonZeros                           = sparse_csr_rowIndices_kae(n_rowIndices)-1

            call mkl_dcsradd('N',0,4,nMatSize,nMatSize, &
                                            sparse_smatrix%csr_values,sparse_smatrix%csr_columns,sparse_smatrix%csr_rowIndices, &
                                          -0.5d0*LF_Aero*LF_Struct,sparse_csr_values_kae(1:nNonZeros),sparse_csr_columns_kae(1:nNonZeros),temp_csr_rowIndices, &
                                            sparse_smatrix_csr_values,sparse_smatrix_csr_columns,sparse_smatrix_csr_rowIndices,nMatSize*nMatSize,info)

          else  !< use only structural sparse matrix
            sparse_smatrix_csr_rowIndices = sparse_smatrix%csr_rowIndices
            sparse_smatrix_csr_columns(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1) = sparse_smatrix%csr_columns
            sparse_smatrix_csr_values(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1)  = sparse_smatrix%csr_values
          end if
          
          !< determine the direction of current aerodynamic force vector and project the aerodynamic force vector at the beginning of the structural step
          if (inos .gt. 1) then 
            vec_norm_w = norm2(the_model_fsi%fae)
            if (vec_norm_w .ne. 0.0d0) then
              unit_fae = the_model_fsi%fae/vec_norm_w
              fae_n_proj = dot_product(fae_n,unit_fae)*unit_fae
            end if
          end if
          
          !< vector of generalized aerodynamic forces for current time step - TO DO - UPDATE THE LINEARIZATION FOR THIS LOAD STEPPING
          if (.not. boolean_dfae) then
            fae = fae_n_proj + (the_model_fsi%fae - fae_n_proj)*LF_Struct
          else !< just a simplification here, it has to be addjusted
            fae = fae_n + (the_model_fsi%fae - fae_n)*LF_Struct
          end if
            
          !< in case of weak coupling, switch off boolean_update_aero
          if (.not. the_model_fsi%boolean_strongfsi) boolean_update_aero = .FALSE.
          
        else !< if fsi-coupling "weak" is activated
          sparse_smatrix_csr_rowIndices = sparse_smatrix%csr_rowIndices
          sparse_smatrix_csr_columns(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1) = sparse_smatrix%csr_columns
          sparse_smatrix_csr_values(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1)  = sparse_smatrix%csr_values
        end if

        !< construction of the residual vector
        rhs_unbalance_force(indicesq) = the_model_structure%fqint + the_model_structure%fqdyn + the_model_structure%fqL - &
                                        the_model_structure%fqgrav - the_model_structure%fqext - &
                                        LF_Aero*fae
                                        
        rhs_unbalance_force(indicesv) = the_model_structure%fv
        rhs_unbalance_force(indicesg) = the_model_structure%g

        !< solving linear system of equation
        nNonZeros = sparse_smatrix_csr_rowIndices(nMatSize+1)-1
        if (boolean_sparse) then
          error = dss_create(handle, MKL_DSS_DEFAULTS)
          error = dss_define_structure(handle, MKL_DSS_NON_SYMMETRIC, sparse_smatrix_csr_rowIndices, nMatSize, nMatSize, sparse_smatrix_csr_columns(1:nNonZeros), nNonZeros )
          error = dss_reorder(handle, MKL_DSS_METIS_OPENMP_ORDER, perm)
        end if
        call solver(nMatSize, handle, sparse_smatrix_csr_values(1:nNonZeros), sparse_smatrix_csr_columns(1:nNonZeros), sparse_smatrix_csr_rowIndices, rhs_unbalance_force, 1, deltaxvector, boolean_sparse, statOut, boolean_determinant)

        !< print determinat
        if (boolean_determinant) then
          inquire(unit=UnIn_det, opened=boolean_itsopen)
          if (.not. boolean_itsopen) then
            call OpenFOutFileToRead(UnIn_det, 'determinant.dres', io_error)
          end if
          write(UnIn_det,'(1000000f30.15)') statOut(1), statOut(2)
        end if

        if (boolean_sparse) then
          error = dss_delete(handle, MKL_DSS_MSG_LVL_WARNING+MKL_DSS_TERM_LVL_ERROR)
        end if

        if (the_model_aero%time .eq. 0.0d0) then
          deltaxvector = 0.0d0
          residuum     = 0.0d0
          exit
        end if

        residuum = dsqrt(dot_product(deltaxvector, deltaxvector))/dsqrt(dot_product(q1, q1)+dot_product(v1, v1)+dot_product(lambda, lambda))
        print '(8X,A9,1X,E10.2E2,1X,A17,1X,E10.2E3)', 'residuum:', residuum, 'convergence work:',abs(dot_product(rhs_unbalance_force, deltaxvector))
        
        if (isnan(residuum)) then
          residuum = 1.0d0
          the_model_structure%boolean_abort = .True.
          exit
        endif

        !< additive actualization of generalized displacements and Lagrange multiplier
        q1     = q1     + deltaxvector(indicesq)
        v1     = v1     + deltaxvector(indicesv)
        lambda = lambda + deltaxvector(indicesg)

        ! check simulation and iteration time
        call system_clock(delta_sim_clock_end_structure)
        delta_sim_time_structure = delta_sim_time_structure + (real(delta_sim_clock_end_structure -delta_sim_clock_start_structure))/real(sim_clock_rate)
      end do
!< ------------- END NEWTON ITERATION -----------------------------------------------------------------------------------------------------    
      !< Abort criteria
      if (residuum > settings%tolerance) then
        print*, '... solution failed to converge at time', time
        if (.not. the_model_fsi%boolean_strongfsi) then
          the_model_structure%boolean_abort = .TRUE.
          exit
        elseif (the_model_fsi%boolean_strongfsi) then
          if (boolean_dfae) then
            the_model_structure%boolean_abort = .TRUE.
            exit
          else
            print*, '... restart load step with linearized aerodynamic forces'
            boolean_dfae = .TRUE.
            istep = istep - 1
            the_model_structure%boolean_abort = .False.
          end if
        end if
      else
        !< Actualization of the old soltion
        q0_test = q0
        q0      = q1
        v0      = v1
        lambda0 = lambda
        print '(5X,A22,1X,F10.5,1X,A11,1X,I2,1X,A9,1X,E10.2E3,1X,A3,F10.5)','converged solution at:',time,'iterations:',iteration,'residuum:',residuum,'t:',delta_sim_time_structure

      end if

    end do
!< ------------- END STRUCTURAL TIME STEPS -----------------------------------------------------------------------------------------------------
     if (the_model_structure%boolean_abort) then
      exit
    end if

    !< update aero grid to current solution
    call the_model_fsi%actualize(the_model_aero, q1, the_model_structure%tncoordinates, v1) ! assign structure's geometry to aerogrid's geometry
    call the_model_aero%actualize(.TRUE., .FALSE.) ! actualize geometry
    call the_model_aero%wake_update_first_row()

    !< update wake data for current aero time step
    the_model_aero%airflow%boolean_actualize = .TRUE. ! this has to be set to true, because the_model_aero%wake_velocity uses free-flow field at wake nodes
    call the_model_aero%wake_velocity(the_model_aero%cutoff)
    the_model_aero%oldcirculations = the_model_aero%circulations

    #if Controller_ADDON
        !< update controller data
        if (the_model_controller%boolean_abort .eqv. .FALSE.) call the_model_controller%update_servo(the_model_aero, the_model_structure, q1, v0_aero, v1)
    #endif

    !< calculate invariants of current time step
    call the_model_structure%modelinvariants(q1, v1)

    !< write results to file. Hint: circulation and loads are determined at time step n+1/2
    call the_model_aero%output_write()
    call the_model_structure%output_write(the_model_aero%tcounter,the_model_aero%time, the_model_aero%time, the_model_aero%tcounter, settings%i_simutype, iteration, q1, v1, lambda, the_model_structure%stress)
    #if Controller_ADDON
        if (the_model_controller%boolean_abort .eqv. .FALSE.) call the_model_controller%output_write()
    #endif

    !< tcounter
    the_model_aero%tcounter = the_model_aero%tcounter+1
    
    !< update fae_n
    fae_n = fae
    
    call system_clock(delta_sim_clock_end_aero)
    delta_sim_time_aero = (real(delta_sim_clock_end_aero - delta_sim_clock_start_aero))/real(sim_clock_rate)
    print '(2X,A31,1X,F10.5,1X,A3,F10.5)','... Aero converged solution at:',the_model_aero%time,'t:',delta_sim_time_aero

  end do 
!< ------------- END AERODYNAMIC TIME STEPS -------------------------------------------------------------------------------------------------------------
  
  if (the_model_structure%boolean_abort) then
    print '(5X,A23,1X,F10.5,1X,A11,1X,I2,1X,A9,1X,E10.2E3)','no convergence at time:',time,'iterations:',iteration,'residuum:',residuum
    goto 999
  end if

  !< Updating model variables
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

  the_model_structure%fqext_t      = the_model_structure%fqext
  the_model_structure%fqext0_t     = fqext0
  the_model_structure%fqint_t      = the_model_structure%fqint
  the_model_structure%fqdyn_t      = the_model_structure%fqdyn
  the_model_structure%fqL_t        = the_model_structure%fqL
  the_model_structure%fv_t         = the_model_structure%fv

  if (settings%output_flag == 1) then
    print *, "Output: Writing matrices..."

    call writeIntVectorToFile(indicesq, 9999,'indicesq.dres')
    call writeIntVectorToFile(indicesv, 9999,'indicesv.dres')
    call writeIntVectorToFile(indicesg, 9999,'indicesg.dres')

    call writeIntVectorToFile(sparse_smatrix_csr_columns(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1), 9999,'csr_columns_S.dres')
    call writeIntVectorToFile(sparse_smatrix_csr_rowIndices,9999,'csr_rowIndices_S.dres')
    call writeRealVectorToFile(sparse_smatrix_csr_values(1:sparse_smatrix_csr_rowIndices(nMatSize+1)-1), 9999,'csr_values_S.dres')

	  call writeIntVectorToFile(the_model_structure%sparse_mass%csr_columns(1:the_model_structure%sparse_mass%csr_rowindices(n_rowIndices)-1),   9999,'csr_columns_M.dres')
    call writeIntVectorToFile(the_model_structure%sparse_mass%csr_rowindices,9999,'csr_rowIndices_M.dres')
    call writeRealVectorToFile(the_model_structure%sparse_mass%csr_values(1:the_model_structure%sparse_mass%csr_rowindices(n_rowIndices)-1),9999,'csr_values_M.dres')

    call writeIntVectorToFile(sparse_csr_columns_kae(1:sparse_csr_rowIndices_kae(n_rowIndices)-1),   9999,'csr_columns_Kae.dres')
    call writeIntVectorToFile(sparse_csr_rowIndices_kae,9999,'csr_rowIndices_Kae.dres')
    call writeRealVectorToFile(sparse_csr_values_kae(1:sparse_csr_rowIndices_kae(n_rowIndices)-1),9999,'csr_values_Kae.dres')

    !!< write dG with in structural coordinates
    !! expand sparse matrix
    !call mkl_ddnscsr(convert_job_csr_to_den, 2*the_model_aero%tncoordinates, 2*the_model_structure%tncoordinates, the_model_fsi%T2n, 2*the_model_aero%tncoordinates, &
    !  the_model_fsi%sparse_T2n%csr_values, the_model_fsi%sparse_T2n%csr_columns, the_model_fsi%sparse_T2n%csr_rowIndices, info)
    !
    !dense_dG = 0.0d0
    !dense_dG = matmul(the_model_aero%dG,the_model_fsi%T2n)
    !call writeRealMatrixToFile(dense_dG,9999,'dG_qv.dres')

  end if

999 continue

  deallocate(amatrix)
  deallocate(deltaxvector)
  deallocate(rhs_unbalance_force)
  deallocate(fqext0)

  return
end subroutine fsi_dynamic_solver
