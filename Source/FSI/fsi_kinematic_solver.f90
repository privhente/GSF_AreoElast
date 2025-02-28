subroutine fsi_kinematic_solver(the_model_structure, the_model_aero, the_model_fsi, settings, boolean_output_write)
! importing classes
  use class_model_structure, only: model_structure, step
  use class_model_aero,      only: model_aero
  use class_model_fsi,       only: model_fsi
  use solver_variables
 
  implicit none

  type(model_structure), intent(inout) :: the_model_structure
  type(model_aero),      intent(inout) :: the_model_aero
  type(model_fsi),       intent(inout) :: the_model_fsi
  type(step),            intent(in)    :: settings
  logical, intent(in) :: boolean_output_write
  integer :: io_error, nsteps
  
  q1 = the_model_structure%q_t
  v1 = the_model_structure%qdot_t
  
  lambda = 0.0d0
  the_model_structure%stress = 0.0d0
  
  !< reading input files from structural FEM simulation
  open(unit = 9999, file = the_model_structure%oldsimfilename // '_qfem.dres', status = 'old', iostat = io_error)
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
  open(unit = 9999, file = the_model_structure%oldsimfilename // '_qfem.dres', status = 'old', iostat = io_error)
  if (io_error .ne. 0) then
      print*, 'Warning: error in reading coordinates file...'
      call exit(-1)
  end if
  
  open(unit = 9998, file = the_model_structure%oldsimfilename // '_vfem.dres', status = 'old', iostat = io_error)
  if (io_error .ne. 0) then
      print*, 'Warning: error in reading coordinates file...'
      !call exit(-1)
  end if

  istep = 0
  do while (istep .le. nsteps)
    istep = istep + 1
    
    !< read kinematic from solution file
    read(9999, *, iostat = io_error) q1
    if (io_error .ne. 0) then
      the_model_structure%boolean_abort = .TRUE.
      print*, 'Warning: error in reading coordinates at step', istep
      exit
    end if
    
    v1 = 0.0d0
    read(9998, *, iostat = io_error) v1
    if (io_error .ne. 0) then
      the_model_structure%boolean_abort = .TRUE.
      print*, 'Warning: error in reading velocities at step', istep
      v1 = 0.0d0
      !exit
    end if

    !< transering structral coordinates to nodes of aerodynamic grid
    call the_model_fsi%actualize(the_model_aero, q1, the_model_structure%tncoordinates, v1) ! assign structure's geometry to aerogrid's geometry
    call the_model_aero%actualize(.TRUE., .FALSE.) ! actualize geometry

    !< writing in result file
    if (boolean_output_write .eqv. .TRUE.) then
      the_model_aero%tcounter = istep
      the_model_aero%time = dble(istep)
      call the_model_aero%output_write()
      call the_model_structure%output_write(istep, dble(istep), dble(istep), istep, 0, istep, q1, v1, lambda, the_model_structure%stress)
    end if
     
    print '(5X,A13,1X,I4)','reading step:', istep

  end do
  
  close(unit = 9999)  
  return
end subroutine fsi_kinematic_solver
