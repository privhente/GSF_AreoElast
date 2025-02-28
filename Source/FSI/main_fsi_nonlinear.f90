subroutine main_fsi_nonlinear()
!< Import classes
  use mkl_dss
  use class_model_structure, only: model_structure
  use solver_functions
  use solver_variables
  use class_model_aero, only: model_aero
  use class_model_fsi, only: model_fsi
  #if Controller_ADDON
    use class_controller, only: model_controller
  #endif
  use my_FileIO_functions

implicit none

!< Defining model_structure and mode_aero
  type(model_structure) :: the_model_structure
  type(model_aero) :: the_model_aero
  #if Controller_ADDON
    type(model_controller) :: the_model_controller
  #endif
  type(model_fsi) :: the_model_fsi

  logical :: boolean_abort = .FALSE., boolean_IO_error = .FALSE.
  real(kind = 8)  :: sim_time = 0.00
  integer(kind=8) :: sim_clock_rate, delta_sim_clock_start, delta_sim_clock_end
  integer(kind=8) :: clock_start, clock_end, clock_rate

  call system_clock(clock_start, clock_rate)

  print*, '-----------------------------------------'
  print*, 'Starting DeSiO-FSI.exe...'
  print*, '-----------------------------------------'

!< Reading input files
  call the_model_structure%readinginputs()
  call the_model_aero%readinginputs()
  #if Controller_ADDON
    call the_model_controller%readinginput()
  #endif
  call the_model_fsi%readinginputs(the_model_structure,the_model_aero)
 
!< Initializing FSI
  if (the_model_fsi%boolean_abort) then
    print*, 'Error: simulation abort due to error in fsi...'
    if (the_model_aero%boolean_abort .eqv. .FALSE.) then
      print*, 'Continuing with DeSiO-Aero.exe...'
      call main_aero(the_model_aero)
      call exit(-1)
    else if (the_model_structure%boolean_abort .eqv. .FALSE.) then
      print*, 'Continuing with DeSiO-Structure.exe...'
      call main_structure(the_model_structure)
      call exit(-1)
    else
      print*, 'Error: simulation stopped...'
      call exit(-1)
    end if
  end if

  print*, 'Initializing simulation model for Fluid-Structur-Interaction'
  the_model_structure%boolean_abort = .FALSE.
  the_model_aero%boolean_abort = .FALSE.

  call the_model_structure%constructor()
  if (the_model_structure%boolean_abort) then
    print*, 'Error: simulation abort due to error in structure...'
    call exit(-1)
  end if

  call the_model_aero%constructor()
  if (the_model_aero%boolean_abort) then
    print*, 'Error: simulation abort due to error in aero...'
    call exit(-1)
  end if

  call the_model_fsi%constructor(the_model_structure,the_model_aero)
  if (the_model_fsi%boolean_abort) then
    print*, 'Error: simulation abort due to error in fsi...'
    call exit(-1)
  end if

  if (the_model_structure%boolean_oldsim) then
    print*, 'Continuing structural simulation from previous result files'
    print*, '... loading previous result files'
    call the_model_structure%read_oldsim()
    #if Controller_ADDON
        the_model_controller%boolean_oldsim = the_model_structure%boolean_oldsim
    #endif
 end if

  if (the_model_aero%boolean_oldsim) then
    print*, 'Continuing aerodynamic simulation from previous result files'
    print*, '... loading previous result files'
    call the_model_aero%read_oldsim()
    print*, '... actualizing aerodynamic geometry and ring circulations'
    if (the_model_aero%boolean_oldsim) then
      call the_model_aero%actualize(.TRUE., .TRUE.)
    else
      the_model_structure%boolean_oldsim = .FALSE.
      #if Controller_ADDON
            the_model_controller%boolean_oldsim = .FALSE.
      #endif
   end if
  end if
  the_model_aero%bool_const_a_matrix = .FALSE.
  the_model_aero%airflow%boolean_actualize = .TRUE.

  call the_model_aero%output_open()
  the_model_structure%boolean_write_init_state = .FALSE.
  call the_model_structure%output_write_model()
  call the_model_structure%output_open()

  #if Controller_ADDON
    if (the_model_controller%boolean_abort .eqv. .FALSE.) call the_model_controller%constructor(the_model_fsi%output_filename, the_model_structure%oldsimfilename, the_model_aero%deltat)
    if (the_model_controller%boolean_abort .eqv. .FALSE.) call the_model_controller%output_open(the_model_fsi%output_filename, the_model_structure%oldsimfilename)
  #endif

  print*, '-------------------------------------------------------'
  print*, 'Running nonlinear FSI...'
  print*, '-------------------------------------------------------'
  print '(1x,A18,1X,F10.3)', 'total time      = ', the_model_aero%totalt
  print '(1x,A18,1X,F10.3)', 'delta time      = ', the_model_aero%deltat
  print '(1x,A18,1X,F10.3)', 'cut off         = ', the_model_aero%cutoff
  print '(1x,A18,1X,A10)',    'simulation type = ', the_model_structure%settings(1)%simutype
  print '(1x,A18,1X,A10)',    'fsi             = ', the_model_fsi%str_fsi_type
  print '(1x,A18,1X,I1)',    'linear. aero    = ', the_model_fsi%flag_dfae
  print*, '-------------------------------------------------------'

  call solver_initialization(the_model_structure)
  the_model_structure%deltat = the_model_structure%settings(1)%deltat

  print*, the_model_structure%boolean_oldsim

  if (the_model_structure%boolean_oldsim .eqv. .FALSE.) then
    !< reading initial state for structure
    print*, 'simulation starts from initial state!'
    !print*, 'reading initial-state from file: ', the_model_structure%oldsimfilename
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

  select case (the_model_structure%settings(1)%simutype)
  case('dynamic')
    call the_model_structure%localSimuSetting(the_model_structure%settings(1)%deltat, 0,.TRUE.,.TRUE.)
    #if Controller_ADDON
        call fsi_dynamic_solver(the_model_structure, the_model_aero, the_model_controller, the_model_fsi, the_model_structure%settings(1),.TRUE.)
    #else
        call fsi_dynamic_solver(the_model_structure, the_model_aero, the_model_fsi, the_model_structure%settings(1),.TRUE.)
    #endif
  case('static')
    call the_model_structure%localSimuSetting(the_model_structure%settings(1)%deltat, 1,.TRUE.,.TRUE.)
    call fsi_static_solver(the_model_structure, the_model_aero, the_model_fsi, the_model_structure%settings(1),.TRUE.)
  case('kinematic')
    call fsi_kinematic_solver(the_model_structure, the_model_aero, the_model_fsi, the_model_structure%settings(1),.TRUE.)
  end select

  !< writing result in check.log file
  allocate(temp_check(1))
  if (the_model_structure%boolean_abort) then
    temp_check = 0.0d0
    call writeRealVectorToFile(temp_check,9999,'check.log')
  else
    temp_check = dsqrt(dot_product(the_model_structure%q_t,the_model_structure%q_t))
    call writeRealVectorToFile(temp_check,9999,'check.log')
  end if
  deallocate(temp_check)

  print*, 'Closing output files...'
  call the_model_aero%output_close()
  call the_model_structure%output_close()
  #if Controller_ADDON
    if (the_model_controller%boolean_abort .eqv. .FALSE.) then
      call the_model_controller%output_close()
    end if
  #endif

  call system_clock(clock_end)
  sim_time = real(clock_end-clock_start) / real(clock_rate)

  if (the_model_structure%boolean_abort .eqv. .FALSE.) then
    print*, ' Finsish DeSiO-FSI...'
    print*, ' Simulation finished'
    print*, ' Wall time = ', sim_time, ' seconds. '
#if B05_ADDON
    flush(6)
    call write_solver_time
#endif
  else
    print*, 'Error: simulation crashed!'
  end if

  return
end subroutine main_fsi_nonlinear
