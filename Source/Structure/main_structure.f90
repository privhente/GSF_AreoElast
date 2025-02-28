! DeSiO main program for structural analysis
subroutine main_structure(the_model_structure)
  use class_model_structure, only: model_structure
  use solver_functions
  use solver_variables
  
  implicit none 
  type(model_structure) :: the_model_structure
  
  real :: start, finish
  real :: clock_rate
  integer(kind=8) :: clock_start, clock_end
  integer :: Simu_Type_Dynamic, Simu_Type_Static, Simu_Type_Buckling, Simu_Type_Modal, Simu_Type_Static_arc, Simu_Type_Invariants, i_simu
  logical :: boolean_IO_error = .FALSE.
  
  Simu_Type_Dynamic    = 0
  Simu_Type_Static     = 1
  Simu_Type_Buckling   = 2
  Simu_Type_Modal      = 3
  Simu_Type_Static_arc = 4
  Simu_Type_Invariants = 5
  
  call system_clock(clock_start, clock_rate)
  call cpu_time(start)

  print*, '-----------------------------------------'
  print*, 'Starting DeSiO-Structure...'
  print*, '-----------------------------------------'

  if (the_model_structure%boolean_model_read .eqv. .FALSE.) then
    print*, 'Reading input files...'
    call the_model_structure%readinginputs()
  end if
  if (the_model_structure%boolean_abort .eqv. .TRUE.) then
    print*, 'Error: simulation aborts in DeSiO-Structure, due to error in or missing files...'
    call exit(-1)
  end if

!< Call constructor, initialization of structural model  
  print*, 'Initializing simulation model...'
  call the_model_structure%constructor()
  call solver_initialization(the_model_structure)

  if (the_model_structure%boolean_oldsim) then
    print*, 'Continuing structure simulation from previous result files'
    print*, '... loading previous result files'
    call the_model_structure%read_oldsim()
  end if
  
  if (the_model_structure%boolean_oldsim .eqv. .FALSE.) then
    !< reading initial state for structure
    print*, 'simulation starts from initial state!'
    print*, 'reading initial-state from file'
    q1     = the_model_structure%q_0
    v1     = the_model_structure%qdot_0
    lambda = 0.0d0
    call readprecalcfiles(q1, v1, lambda, the_model_structure%oldsimfilename, boolean_IO_error)
    the_model_structure%q_t      = q1
    the_model_structure%qdot_t   = v1
    the_model_structure%lambda_t = lambda
    if (boolean_IO_error) then
      print*, '... reading initial-state from file failed due to missing files!'
    else
      print*, '... succesfully read initial-state!'
    end if
  endif
  
  call the_model_structure%output_write_model()
  
  ! Loop over all simulation steps
  do i_simu = 1, size(the_model_structure%settings)
    print *, "running simulation with timestep: ", the_model_structure%settings(i_simu)%deltat
    print *, "Simulation type: ", the_model_structure%settings(i_simu)%simutype

    select case (the_model_structure%settings(i_simu)%simutype)
      
    case ('dynamic')
        if (.not. the_model_structure%boolean_outputfiles_open) call the_model_structure%output_open()
        call the_model_structure%localSimuSetting(the_model_structure%settings(i_simu)%deltat, Simu_Type_Dynamic,.TRUE.,.TRUE.)
        the_model_structure%settings(i_simu)%i_simutype = Simu_Type_Dynamic
        call dynamic_solver(the_model_structure, the_model_structure%settings(i_simu),.TRUE.)
        
    case ('static')
        if (.not. the_model_structure%boolean_outputfiles_open) call the_model_structure%output_open()
        call the_model_structure%localSimuSetting(the_model_structure%settings(i_simu)%deltat, Simu_Type_Static,.TRUE.,.TRUE.)
        the_model_structure%settings(i_simu)%i_simutype = Simu_Type_static
        call static_solver(the_model_structure, the_model_structure%settings(i_simu),.TRUE.)

    case ('static_arc')
        if (.not. the_model_structure%boolean_outputfiles_open) call the_model_structure%output_open()
        call the_model_structure%localSimuSetting(the_model_structure%settings(i_simu)%deltat, Simu_Type_Static_arc,.TRUE.,.TRUE.)
        the_model_structure%settings(i_simu)%i_simutype = Simu_Type_Static_arc
        call static_solver_arc(the_model_structure, the_model_structure%settings(i_simu),.TRUE.)
        
    case ('buckling')
        if (.not. the_model_structure%boolean_outputfiles_open) call the_model_structure%output_open()
        the_model_structure%settings(i_simu)%gravityflag = 0
        call the_model_structure%localSimuSetting(1.0d0, Simu_Type_Buckling,.TRUE.,.TRUE.)
        the_model_structure%settings(i_simu)%i_simutype = Simu_Type_Buckling
        call buckling_solver(the_model_structure, the_model_structure%settings(i_simu),.TRUE.)
        
    case ('modal')
        if (.not. the_model_structure%boolean_outputfiles_open) call the_model_structure%output_open()
        the_model_structure%settings(i_simu)%gravityflag = 0
        call the_model_structure%localSimuSetting(1.0d0, Simu_Type_Modal,.TRUE.,.TRUE.)
        the_model_structure%settings(i_simu)%i_simutype = Simu_Type_modal
        call modal_solver(the_model_structure, the_model_structure%settings(i_simu),.TRUE.)

    case ('invariants')
        if (the_model_structure%boolean_outputfiles_open) then
          call the_model_structure%output_close(6)
          call the_model_structure%output_close(7)
          call the_model_structure%output_close(8)
        end if
        the_model_structure%boolean_write_init_state = .FALSE.
        call the_model_structure%output_open(6)
        call the_model_structure%output_open(7)
        call the_model_structure%output_open(8)
        !call the_model_structure%localSimuSetting(1.0d0, Simu_Type_Invariants,.TRUE.,.TRUE.)
        call invariants(the_model_structure, the_model_structure%settings(i_simu), .TRUE.)
    end select

    if (the_model_structure%boolean_abort) then
      print*,'Simulation crashed !!'
      exit
    end if
    
  end do
  
! Closing result files  
  print*, 'Closing output files...'
  call the_model_structure%output_close()
  
  call system_clock(clock_end)
  call cpu_time(finish)
  
  print*, ' Simulation finished'
  print*, ' CPU  time = ', finish-start, ' seconds. '
  print*, ' Wall time = ', real(clock_end-clock_start) / clock_rate, ' seconds. '    
      
  if (the_model_structure%boolean_abort) then
    stop
  end if
          
end subroutine main_structure