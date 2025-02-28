subroutine main_aero(the_model_aero)
!< Import classes
  use mkl_dss
  use class_model_aero, only: model_aero
  use my_FileIO_functions, only : writeRealVectorToFile, writeRealMatrixToFile
  
implicit none

!< Defining model_structure and mode_aero  
  type(model_aero):: the_model_aero
  real(kind = 8)  :: sim_time = 0.00
  integer(kind=8) :: sim_clock_rate, delta_sim_clock_start, delta_sim_clock_end
  integer(kind=8) :: clock_start, clock_end, clock_rate
  real(kind = 8), allocatable :: temp_check(:)
  
  call system_clock(clock_start, clock_rate)
  print*, '-----------------------------------------'
  print*, 'Starting DeSiO-Aero...'
  print*, '-----------------------------------------'
  
  if (the_model_aero%boolean_model_read .eqv. .FALSE. ) then
    print*, 'Reading input files...'
    call the_model_aero%readinginputs()
  end if
  if (the_model_aero%boolean_abort .eqv. .TRUE.) then
    print*, 'Error: simulation abort in DeSiO-Aero, due to missing aerodynamic files...'
    stop
  end if

!< Call constructor, initialization of aerodynamic model
  print*, 'Initializing simulation model...'
  call the_model_aero%constructor()
  
!< Initializing FSI
  print*, '-------------------------------------------------------'
  print*, 'Running calculation with explicit time integration...'
  print*, '-------------------------------------------------------'
  print*, 'total time = ', the_model_aero%totalt
  print*, 'delta time = ', the_model_aero%deltat
  print*, 'cut off    = ', the_model_aero%cutoff
  print*, '-------------------------------------------------------'
  
  call the_model_aero%read_kinematic()
  
  if (the_model_aero%boolean_oldsim) then
    print*, 'continuing simulation from previous result files'
    print*, 'reading previous result files'
    call the_model_aero%read_oldsim()
    if (the_model_aero%boolean_oldsim) then
      print*, 'actualize geometry and ringcirculations'
      call the_model_aero%actualize(.TRUE., .TRUE.)
    end if
  end if
  
  !< set airflow flag to TRUE, so that external flow field is updated in each aero time step
  the_model_aero%airflow%boolean_actualize = .TRUE.
  the_model_aero%boolean_unsteady_term = .TRUE.
  
  !< open results files for output
  call the_model_aero%output_open()

  if (the_model_aero%strSimuType .ne. 'num_sensitivity') then
    
    print*, 'Starting aerodynamical simulation'
    !< Loop over aerodynamic time step
    do while (the_model_aero%tcounter+1 .le. the_model_aero%nsteps)
      call system_clock(delta_sim_clock_start, sim_clock_rate)
        
      the_model_aero%time = the_model_aero%tcounter*the_model_aero%deltat
      print '(2X,A10,1X,I4,1X,A5,1X,F10.3,1X,A3,1X,F10.3)', 'Aero step:', the_model_aero%tcounter, 'time:', the_model_aero%time, 'dt:', the_model_aero%delta_sim_time
      
      !< performing aerodynamic step
      call the_model_aero%solver()
      call the_model_aero%output_write()
      
      the_model_aero%tcounter = the_model_aero%tcounter+1
    
      !< for output window
      call system_clock(delta_sim_clock_end)
      the_model_aero%delta_sim_time = (real(delta_sim_clock_end -delta_sim_clock_start))/real(sim_clock_rate)
      sim_time  = sim_time + the_model_aero%delta_sim_time
    end do
      
    !< Checking linearized equation of last simulation step
    select case (the_model_aero%strSimuType)
    case ('num_diff_DG')
      print*, 'Checking linearized non-penetration vs. numerical differentation'
      call num_diff_DG(the_model_aero)
    case ('num_diff_q')
      print*, 'Checking dfae,q vs. numerical differentation'
      call num_diff_q(the_model_aero)
    case ('num_diff_v')
      print*, 'Checking dfae,v vs. numerical differentation'
      call num_diff_v(the_model_aero)
    end select    
  
  elseif (the_model_aero%strSimuType .eq. 'num_sensitivity') then

    print*, 'Starting numerical force calculation with frozen wake'
    call system_clock(delta_sim_clock_start, sim_clock_rate)
      
    !< read only wake from file
    call the_model_aero%read_oldwake()
    
    !< for numerical calculation of one step with frozen wake - in steady- state
    the_model_aero%time = the_model_aero%tcounter*the_model_aero%deltat
    call the_model_aero%solver_num_sensitivity()
      
    print '(2X,A10,1X,I4,1X,A5,1X,F10.3,1X,A3,1X,F10.3)', 'Aero step:', the_model_aero%tcounter, 'time:', the_model_aero%time, 'dt:', the_model_aero%delta_sim_time
    call the_model_aero%output_write()
      
    !< for output window
    call system_clock(delta_sim_clock_end)
    the_model_aero%delta_sim_time = (real(delta_sim_clock_end -delta_sim_clock_start))/real(sim_clock_rate)
    sim_time  = sim_time + the_model_aero%delta_sim_time
    
  end if
  
  !< writing result in check.log file
  allocate(temp_check(1))
  if (the_model_aero%boolean_abort) then
    temp_check = 0.0d0
    call writeRealVectorToFile(temp_check,9999,'check.log')
  else
    temp_check = dsqrt(dot_product(the_model_aero%circulations,the_model_aero%circulations))
    call writeRealVectorToFile(temp_check,9999,'check.log')
  end if
  deallocate(temp_check)
    
  !< Closing ouptut files  
  print*, 'Closing output files...'
  call the_model_aero%output_close()
  
  call system_clock(clock_end)
  sim_time = real(clock_end-clock_start) / real(clock_rate)
    
  print*, ' Finsish DeSiO-Aero...'
  print*, ' Simulation finished'
  print*, ' Wall time = ', sim_time, ' seconds. ' 
  
  return
  end subroutine main_aero
