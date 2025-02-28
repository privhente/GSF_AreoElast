module class_controller
    use class_model_structure,  only: model_structure
    use class_model_aero,       only: model_aero
    use class_inflow_aero,      only: air_flow
    use my_math_structure,      only: cross
    use my_constants_structure, only: pi
    use my_FileIO_functions,    only: GetNewUnit
    implicit none
    
    type :: model_controller
        character(len = 7)                        :: softwaretype 
        character(len=:), allocatable             :: infile
        real(kind = 8), dimension(3)              :: v1_tt, v2_tt, n_shaft, azimuth  
        real(kind = 8)                            :: init_pitch, swap(160)
        real(kind = 8), allocatable, dimension(:) :: pitch_t, pitch_t0, pitch_0, M_IPB, M_OPB, assign_pitch_t
        real(kind = 8)                            :: momega_hub, omega_hub, thrust, v_wind, mgenerator
        real(kind = 8)                            :: time, deltat, time_init
        real(kind = 8)                            :: p_el 
        real(kind = 8)                            :: shaft_pwr, m_rotor_aero, del_phi(3), hub_aux_1(12), hub_aux_2(12)
        integer, allocatable, dimension(:)        :: revolutej2_ID, rotationl_ID, blade_nbr
        integer                                   :: number_blades, length_infile
        integer                                   :: UnIn_servofile
        logical                                   :: boolean_abort = .FALSE., boolean_oldsim = .FALSE.
        real(kind = 8), dimension(3)              :: M_OPB_2
        
    contains
    
        procedure :: constructor
        procedure :: assign_servo_to_model
        procedure :: update_servo
        procedure :: readinginput
        procedure :: output_open
        procedure :: output_write
        procedure :: output_close
        procedure :: read_oldsim
    end type model_controller
        
  contains
  
  !< subroutine to initialize and allocate some controller parameter and arrays
  subroutine constructor(this, output_filename, oldsimfilename, deltat)
    implicit none

    class(model_controller), intent(inout) :: this
    character(:), allocatable :: output_filename, oldsimfilename
    integer :: i
    integer :: UnIn_servofile, io_error, nrows
    real(kind = 8) :: deltat
    logical :: boolean_opened
    
    allocate(this%pitch_t(this%number_blades))
    allocate(this%pitch_0(this%number_blades))
    allocate(this%pitch_t0(this%number_blades))
    allocate(this%assign_pitch_t(this%number_blades))
    allocate(this%M_IPB(this%number_blades))
    allocate(this%M_OPB(this%number_blades))

    this%time       = 0.0d0
    this%deltat     = deltat
    this%time_init  = this%deltat
    this%p_el       = 0.0d0
    this%pitch_0    = this%init_pitch*pi/180        !< initial pitch angle given in Input
    this%pitch_t0   = this%pitch_0                  !< variable to store the pitch angle of t-1, for t=0: initial pitch
    this%pitch_t    = this%pitch_0                  !< pitch angle given by controller to be set for next time step
    this%mgenerator = 0.0d0
    this%omega_hub  = 0.0d0
    this%momega_hub = 0.0d0
    this%thrust     = 0.0d0 
    this%v_wind     = 0.0d0
    this%M_OPB      = 0.0d0
    this%M_IPB      = 0.0d0
    
    this%swap       = 0.0d0
    if (this%boolean_oldsim)  then
      print*, 'Continuing servo data from previous result files'
      print*, '... loading previous result files'
      call this%read_oldsim(output_filename, oldsimfilename)
    end if
      
  return
  end subroutine

!< subroutine to read results of previous converged simulation, for on Rosco-controller (libdiscon.dll)
subroutine read_oldsim(this, output_filename, oldsimfilename)
    implicit none

    class(model_controller),  intent(inout) :: this
    character(:), allocatable :: output_filename, oldsimfilename
    integer :: i
    integer :: UnIn_servofile, io_error, nrows
    real(kind = 8) :: deltat
    logical :: boolean_opened
    
1000 format(1000000f30.15)
1100 format(1000000E27.18E3)

    ! read old simulation file
    call GetNewUnit (UnIn_servofile)
    !< reading result files from previous simulation
    call GetNewUnit (UnIn_servofile)
    open(unit = UnIn_servofile, file = oldsimfilename // '_servo.dres', status = 'old', iostat = io_error)
    if (io_error .ne. 0) then
      this%boolean_abort  = .TRUE.
      this%boolean_oldsim = .FALSE.
      print*, '... warning: could not open ', oldsimfilename // '_servo.dres'
    else
      !< determine number of rows in result files
      nrows = 0
      read(UnIn_servofile, *, iostat = io_error)
      do 
        if (io_error .ne. 0) exit
        read(UnIn_servofile, *, iostat = io_error)
        nrows = nrows + 1
      end do
      close(unit = UnIn_servofile)
    end if
      
    call GetNewUnit (UnIn_servofile)
    open(unit = UnIn_servofile, file = oldsimfilename // '_servo.dres', status = 'old', iostat = io_error)
      do i = 1, nrows-2
        read(UnIn_servofile, 1100, iostat = io_error)
      end do
      read(UnIn_servofile, 1100, iostat = io_error) this%time, this%pitch_t0(:), this%mgenerator, this%omega_hub,  this%momega_hub, this%thrust, this%v_wind
      read(UnIn_servofile, 1100, iostat = io_error) this%time, this%pitch_t(:),  this%mgenerator, this%omega_hub,  this%momega_hub, this%thrust, this%v_wind
      if (io_error .ne. 0) then
        this%boolean_abort  = .TRUE.
        this%boolean_oldsim = .FALSE.
        print*, '... warning: error in reading ', oldsimfilename // '_servo.dres'
      else
        print*, '... servo file ', oldsimfilename // '_servo.dres', ' read correctly'
      end if
    close(unit = UnIn_servofile)
    this%time_init = this%time + this%deltat
    
    return
end subroutine read_oldsim

!< subroutine to assign controller data to the structural model
subroutine assign_servo_to_model(this, this_structure, q_2, q_1, time, deltat)
    implicit none

    class(model_controller),  intent(inout) :: this
    class(model_structure),   intent(inout) :: this_structure
    real(kind = 8),           intent(inout) :: q_2(:)
    real(kind = 8),           intent(inout) :: q_1(:)
    
    real(kind = 8)                :: time, deltat, controller_initialization_time
    real(kind = 8)                :: hub_q(12), nac_q(12), bla_q(12), phi(3,2), dir_rj2(3)
    real(kind = 8)                :: T(3,3)
    integer                       :: i, nn_bla, nn_nac, nn_hub
    logical                       :: controller_exists
    
    !< update controller time
    this%time   = time
    this%deltat = deltat
    
    !< rotation axis rotor
    nn_nac = this_structure%constraints12(this%revolutej2_ID(1)-1)%c%nodes(1)  ! revolutej2_ID entspr. rel_rot_locID
    nn_hub = this_structure%constraints12(this%revolutej2_ID(1)-1)%c%nodes(2)
    
    nac_q = q_2(this_structure%nodes12(nn_nac)%coordinates)
    hub_q = q_2(this_structure%nodes12(nn_hub)%coordinates)
    phi   = this_structure%constraints12(this%revolutej2_ID(1)-1)%c%constraint%phi    
    
    this%n_shaft      = 0.0d0
    this%n_shaft      = (hub_q(1:3)+phi(1, 2)*hub_q(4:6)+phi(2, 2)*hub_q(7:9)+phi(3, 2)*hub_q(10:12)) & 
                        - (nac_q(1:3)+phi(1, 1)*nac_q(4:6)+phi(2, 1)*nac_q(7:9)+phi(3, 1)*nac_q(10:12))
    this%n_shaft(1:3) = -this%n_shaft/norm2(this%n_shaft) !< Inverting the rotation axis: hub -> nacelle: Positive rotor rotation 
    
    !< set generator moment: momega_hub determined from last time step - mgenerator determined from controller
    this%mgenerator = this%swap(47)                                                                       !< demanded generator torque that slows down the rotation of the shaft
        !< n_shaft: global -> local
    this_structure%loads12(1)%material(4) = dot_product(this%n_shaft(1:3),hub_q(4:6))
    this_structure%loads12(1)%material(5) = dot_product(this%n_shaft(1:3),hub_q(7:9))
    this_structure%loads12(1)%material(6) = dot_product(this%n_shaft(1:3),hub_q(10:12))
    !this_structure%loadamplitudes12(1)%intensity = 1.0d0 * this%mgenerator * sign(1.0d0,this%momega_hub) !< generator moment always against shaft rotation direction
    this_structure%loadamplitudes12(1)%intensity = -abs(this%mgenerator)
    
    !< generator moment acting in opposite direction on the nacelle
    this_structure%loads12(2)%material(4) = dot_product(this%n_shaft(1:3),nac_q(4:6))
    this_structure%loads12(2)%material(5) = dot_product(this%n_shaft(1:3),nac_q(7:9))
    this_structure%loads12(2)%material(6) = dot_product(this%n_shaft(1:3),nac_q(10:12))
    !this_structure%loadamplitudes12(2)%intensity = -1.0d0 * this%mgenerator * sign(1.0d0,this%momega_hub)  !< sign(a,b) returns value a with sign of b
    this_structure%loadamplitudes12(2)%intensity = abs(this%mgenerator)
    
    !< adapting the blade positions to the new pitch angle
    do i = 1, this%number_blades
      if (this%swap(1) .ne. 0) then 
        this%pitch_t0(i) = this%pitch_t(i)              !< store the last pitch angle
        this%pitch_t(i)  = this%swap(42+(i-1))          !< use the new pitch angle from controller output
      end if  
    
    !< set the controller pitch to the relative_rotation_local (pitch) constraint 
      this%assign_pitch_t(i) = -this%pitch_t(i)   
      this_structure%constraints12(this%revolutej2_ID(i))%c%constraint%dir(1) = this%assign_pitch_t(i)+this%pitch_0(i)

    end do
    
  return
end subroutine assign_servo_to_model

!< subroutine to update control parameter depending on structural response of previous controller input
subroutine update_servo(this, this_aero, this_structure, q_2, v_1, v_2) 
    use, intrinsic :: iso_c_binding, only : c_float, c_int, c_char, c_null_char
#if Controller_ADDON

#if System_Windows
    use DISCON_m, only : SetupDISCON, DISCON
#else
    use DISCON_m_linux, only : DISCON
#endif

#endif    
    implicit none
    
    class(model_controller), intent(inout) :: this
    class(model_structure), intent(inout)  :: this_structure
    class(model_aero),      intent(inout)  :: this_aero

    real(kind = 8), intent(inout)  :: v_1(:), v_2(:), q_2(:)
    integer                        :: nn_tto, r1
    integer                        :: nn_bla, nn_nac, nn_hub
    integer                        :: i
    real(kind = 8), dimension(12)  :: hub_q, hub_v, nac_q, bla_q
    real(kind = 8), dimension(3)   :: fae_total, M_b, ax
    real(kind = 8)                 :: phi(3,2), dir_rj2(3), M_root_gl(3),  forces_bla_gl(3), moments_bla_gl(3)
    real(kind = 8)                 :: a_fa, a_ss, a_tt(3)
    real(kind = 8), dimension(3)   :: v1_tt, v2_tt
    real(kind = 8), dimension(3)   :: M_b_gl, M_b_loc, F_b_gl, F_b_loc
    real(kind = 8), dimension(3)   :: ax_hor, ax_vert, ax_blade, forces_bla_loc, moments_bla_loc
    real(kind = 8), dimension(3,3) :: T
    real(kind = 8), dimension(12)  :: bla_q_1, bla_q_2
    real(kind = 8), dimension(3)   :: m_hub_aero, dist, f_ring
    integer                        :: n_surface, n_ring
        
    !< input variables of DISCON
    real(c_float), dimension(160) :: ctrl_avrSWAP           ! 85 m\FCssten aussreichend sein, siehe dbg3 fast
    integer(c_int)                :: ctrl_aviFAIL     
    character(kind=c_char,len=14) :: ctrl_accINFILE
    character(kind=c_char,len=7)  :: ctrl_avcOUTNAME
    character(kind=c_char,len=49) :: ctrl_avcMSG
    
    if (this%time .gt. 0.0d0) then                                                                      ! First simulation step without controller to initialize all variables
      !< positions and velocities nacelle determined using node numbers of revolutejoint, which will be always one before revolutejoint_2
      nn_nac = this_structure%constraints12(this%revolutej2_ID(1)-1)%c%nodes(1)
      nn_hub = this_structure%constraints12(this%revolutej2_ID(1)-1)%c%nodes(2)
    
      nac_q = q_2(this_structure%nodes12(nn_nac)%coordinates)
      hub_q = q_2(this_structure%nodes12(nn_hub)%coordinates)
      hub_v = v_2(this_structure%nodes12(nn_hub)%coordinates)
        
      !< new (normed) rotation axis of rotor
      this%n_shaft(1:3) = 0.0d0
      this%n_shaft(1:3) = -this_structure%constraints12(this%revolutej2_ID(1)-1)%c%constraint%dir_gl(1:3)! see class_constraint_12.f90 - grevolutejoint, normed
            
      !< calculating angular velocity of rotor
      this%omega_hub = dot_product(0.5d0*(cross(hub_q(4:6),hub_v(4:6))+cross(hub_q(7:9),hub_v(7:9))+cross(hub_q(10:12),hub_v(10:12))), this%n_shaft)
        
      !< calculating rotor moment and in-plane and out-of-plane bending moments at blade root
      M_root_gl = 0.0d0
      do i = 1, this%number_blades
        nn_bla  = this_structure%constraints12(this%revolutej2_ID(i))%c%nodes(2)
        bla_q_1   = q_2(this_structure%nodes12(nn_bla)%coordinates)     !< Coordinates blade root
        bla_q_2   = q_2(this_structure%nodes12(nn_bla+1)%coordinates)   !< Coordinates blade root +1
        bla_q     = 0.5d0*(bla_q_2 + bla_q_1)                           !< Gau\DF-point: Where the forces and moments are evaluated, difference 1/2 element lenght
        
        phi     = this_structure%constraints12(this%revolutej2_ID(i))%c%constraint%phi  !< TODO: check if phi has to be modified
        !dir_rj2 = (bla_q(1:3)+phi(1,2)*bla_q(4:6)+phi(2,2)*bla_q(7:9)+phi(3,2)*bla_q(10:12)) - (hub_q(1:3)+phi(1,1)*hub_q(4:6)+phi(2,1)*hub_q(7:9)+phi(3,1)*hub_q(10:12))       
        dir_rj2 = bla_q_1(10:12)                                        !< Rotation axis blade in golbal coordinates
        
        !< M = sum(M) + sum(rxF) --> checked (moment on the hub) 
        forces_bla_gl  = bla_q(4:6)*this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(1) + &
                         bla_q(7:9)*this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(2) + &
                         bla_q(10:12)*this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(3)
        
        moments_bla_gl = bla_q(4:6)*this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(4) + &
                         bla_q(7:9)*this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(5) + &
                         bla_q(10:12)*this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(6)
        
        M_root_gl  = M_root_gl + moments_bla_gl  + cross(dir_rj2,forces_bla_gl) ! global moment vector on the hub
        
                        !< To check the one above: Another way to calculate the global forces and moments
                        !T(1,1:3) = bla_q(4:6) ! Transformationsmatrix, T^T = T^-1, weil orthogonal
                        !T(2,1:3) = bla_q(7:9)
                        !T(3,1:3) = bla_q(10:12)
                        !
                        !forces_bla_loc(1)  = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(1)
                        !forces_bla_loc(2)  = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(2)
                        !forces_bla_loc(3)  = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(3)
                        !moments_bla_loc(1) = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(4)
                        !moments_bla_loc(2) = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(5)
                        !moments_bla_loc(3) = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(6)
                        !
                        !forces_bla_gl(1:3) = matmul(transpose(T),forces_bla_loc)
                        !moments_bla_gl(1:3) = matmul(transpose(T),moments_bla_loc)

        !< Vector of moments on the blade root (local COS) 
        M_b_loc   = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(4:6) 
        F_b_loc   = this_structure%beams(this%blade_nbr(i))%elements(1)%stress_resultants(1:3)  !< local  forces on midpoint of first element
        F_b_gl    = F_b_loc(1)*bla_q(4:6) + F_b_loc(2)*bla_q(7:9) + F_b_loc(3)*bla_q(10:12)     !< global forces on midpoint of  first element 
        M_b_gl    = M_b_loc(1)*bla_q(4:6) + M_b_loc(2)*bla_q(7:9) + M_b_loc(3)*bla_q(10:12) + cross(0.5d0*(bla_q_2(1:3) - bla_q_1(1:3)), F_b_gl)
    
        !< calculation a perpendicular vector to the pitch axis in the rotor plane to determine the OP bending moment for pitched blades
        ax = cross(this%n_shaft(1:3), dir_rj2)
        ax = ax/norm2(ax)
        this%M_OPB(i) = dot_product(M_b_gl,ax) 
      
        !< TODO: Calculation of in plane bending moments correct?
        !< in plane um die shaft Achse oder je um die Pitch-Achse (this%M_IPB(i) = dot_product(M_b_gl, dir_rj2))
        this%M_IPB(i) = dot_product(M_b_gl-this%swap(47), this%n_shaft)
      
      end do
      
      !< rotor torque
      this%n_shaft = this%n_shaft/norm2(this%n_shaft)
      this%momega_hub = dot_product(M_root_gl,this%n_shaft)
      
     !< calculation the rotor moment based on the aerodynamic forces
      m_hub_aero        = 0.0d0
      this%m_rotor_aero = 0.0d0
      
      do n_surface = 1, this%number_blades*2                                               !< loop over all surfaces
        do n_ring   = 1, this_aero%surfaces(n_surface)%vortex_sheet%nrings                 !< loop over all rings of a surface
          dist(1:3) = this_aero%surfaces(n_surface)%vortex_sheet%rings(n_ring)%center(1:3) !< vector to the center of surface ring
          dist(1:3) = dist(1:3) - hub_q(1:3)                                               !< distance between ring center and hub
          
          f_ring     = this_aero%surfaces(n_surface)%vortex_sheet%rings(n_ring)%fae(1:3)   !< force vector on center of surface ring
          m_hub_aero = m_hub_aero + cross(dist(1:3),f_ring(1:3))                           !< resulting moment on the hub
        end do
      end do
      this%m_rotor_aero = dot_product(m_hub_aero, this%n_shaft)                            !< projecting moment to the rotor axis
      
      !< calculating tower-top accelerations - a = (v_n+1 - v_n)/ (dt)
      nn_tto = this_structure%constraints12(this%revolutej2_ID(1)-2)%c%nodes(1)
      !nn_tto = this_structure%constraints12(this%revolutej2_ID(1)-5)%c%nodes(1)
      v1_tt  = v_1(this_structure%nodes12(nn_tto)%coordinates(1:3))
      v2_tt  = v_2(this_structure%nodes12(nn_tto)%coordinates(1:3))
      a_tt   = (v2_tt(1:3) - v1_tt(1:3))/this%deltat                                                                         
      !< changing to tower-top director cos
      a_fa   = dot_product(a_tt,nac_q(4:6))                                                           ! fore-aft acceleration   
      a_ss   = dot_product(a_tt,nac_q(7:9))                                                           ! side-side acceleration 
            
      !< wind speed at hub height
      this%v_wind = norm2(this_aero%airflow%airflowvinfinity(hub_q(1:3), this%time))
    
      !< calculating rotor thrust force
      fae_total = 0.0d0
      do i = 1,this_aero%tnrings
        fae_total = fae_total + this_aero%fae(3*(i-1)+1:3*(i-1)+3)  !global
      end do 
      this%thrust = dot_product(fae_total, this%n_shaft)
      
      !< electrical power output
      this%p_el = this%mgenerator * this%omega_hub
      
      !< measured shaft power: torque * rotation speed
      !this%shaft_pwr = this%momega_hub*this%omega_hub
      this%shaft_pwr = this%m_rotor_aero*this%omega_hub
      
      !< calculation of the azimuth angle
      !< TODO: define the variables globally that are used multiple times
      nac_q = q_2(this_structure%nodes12(nn_nac)%coordinates)
      ax_hor(:) = nac_q(7:9)
      ax_vert(:) = cross(this%n_shaft, ax_hor) 
      ax_vert(:) = ax_vert/norm2(ax_vert)
      do i = 1, this%number_blades
        nn_bla  = this_structure%constraints12(this%revolutej2_ID(i))%c%nodes(2)
        bla_q   = q_2(this_structure%nodes12(nn_bla)%coordinates)
        phi     = this_structure%constraints12(this%revolutej2_ID(i))%c%constraint%phi
        ax_blade = (bla_q(1:3)+phi(1,2)*bla_q(4:6)+phi(2,2)*bla_q(7:9)+phi(3,2)*bla_q(10:12)) - (hub_q(1:3)+phi(1,1)*hub_q(4:6)+phi(2,1)*hub_q(7:9)+phi(3,1)*hub_q(10:12))
        ax_blade = ax_blade/norm2(ax_blade) 
        ax_blade = bla_q(10:12)
        if (ax_vert(1)*ax_blade(2)-ax_vert(2)*ax_blade(1) <= 0) then                    
          this%azimuth(i) = acos(dot_product(ax_vert,ax_blade)/(norm2(ax_vert)*norm2(ax_blade)))
        else
          this%azimuth(i) = 2*pi - acos(dot_product(ax_vert,ax_blade)/(norm2(ax_vert)*norm2(ax_blade))) !< because angle ranges between 0 and 2pi
        end if
      end do
      
      !< input positions according to bladed user manual
      !call ListDisconInfo()
      if (this%time .eq. this%deltat) then 
          this%swap(1) = 0.0d0
      elseif (this%time .gt. this%deltat) then
          this%swap(1) = 1.0d0
      elseif (this%time .eq. this_aero%totalt) then 
          this%swap(1) = -1.0d0    
      end if
                                                     
      this%swap(2)  = this%time                      !< simulation time
      this%swap(3)  = this%deltat                    !< delta t
      this%swap(4)  = this%pitch_t(1)                !< pitch angle of blade 1
      this%swap(10) = 0.0d0                          !< 0: pitch position actuator, 1 = pitch rate actuator
      !this%swap(11) = this%swap(45)                 !< swap(11): Current demanded pitch, swap(45): Demanded pitch angle, collective 
      this%swap(14) = this%shaft_pwr                 !< Measured shaft power
      this%swap(15) = this%p_el                      !< measured electrical power output
      this%swap(20) = this%omega_hub                 !< generator speed (direct drive has transmission factor 1 => generator speed = rotor speed) 
      !this%swap(21) = this%omega_hub                !< rotor speed
      this%swap(22) = this%swap(47)                  !< no generator model: measured generator torque = demanded generator torque || not included in the FAST DBG3 output
      this%swap(23) = this%swap(47)                  !< no generator model: measured generator torque = demanded generator torque
      this%swap(24) = 0.0d0                          !< measured yaw error 
      this%swap(27) = this%v_wind                    !< horizontal wind speed
      this%swap(28) = 0.0d0                          !< pitch controll: 0 = collective, 1 = individual
      do i = 1, this%number_blades   
        this%swap(30+(i-1)) = -this%M_OPB(i)         !< root out of plane bending moment: this%swap(30), this%swap(31), this%swap(32) = M32 
      end do                                         !< changed sign in order to get positive input values as in the fast input
      do i = 2, this%number_blades                   
        this%swap(33+(i-2)) = this%pitch_t(i)       !< blade 2 blade 3 pitch angles this%swap(33), this%swap(34)
      end do                                   
      this%swap(35) = 1.0d0                          !< Generator contactor
      !this%swap(37) = this%thrust                   !< additional in B02 controller, not included in bladed user manual as input value
      this%swap(49) = 10                             !< Maximum no. of characters allowed in the \93MESSAGE\94, DLL case only   
      this%swap(50) = this%length_infile             !< No. of characters in the \93INFILE\94 argument
      this%swap(51) = 5.0d0                          !< No. of characters in the "OUTNAME" argument   
      this%swap(53) = a_fa                           !< fore-aft acceleration
      this%swap(54) = a_ss                           !< side-side acceleration
      this%swap(60) = this%azimuth(1)                !< rotor azimuth, here blade 1, not clear if a certain blade is demanded
      this%swap(61) = this%number_blades             !< number of blades
      this%swap(62) = 300.0d0                        !< Max. number of values which can be returned for logging
      this%swap(63) = 165.0d0                        !< Record number for start of logging output
      this%swap(64) = 12601.0d0                      !< Max. no. of characters which can be returned in \93OUTNAME\94 
      this%swap(66) = 1001.0d0                       !< reserved input variable - value adapted from openFast input, dbg3 file
      do i = 1, this%number_blades                    
        this%swap(69+(i-1)) = this%M_IPB(i)          !< root in plane bending moment: this%swap(69), this%swap(70), this%swap(71) 
      end do
      
      !< Filling the controller function variables
      if (this_aero%time .eq. this%deltat) then
#if Controller_ADDON

#if System_Windows
        call SetupDISCON(this%boolean_abort)         !< has to be loaded only once in the beginning for windows
#endif

#endif
      end if

      if (this%boolean_abort .eqv. .TRUE.) then
        print*, '... error in SetupDiscon: continuing without controller...'
        return
      end if
      
      !< further inputs of the DISCON controller routine
      ctrl_accINFILE   = this%infile
      ctrl_aviFAIL     = 0_c_float
      ctrl_avcOUTNAME  = c_char_'servo'// c_null_char
      ctrl_avrSWAP     = this%swap                  !< converting input values to C++ float values
      
      !< Calling the actual controller function from the DLL
#if Controller_ADDON      
      call DISCON(ctrl_avrSWAP, ctrl_aviFAIL, ctrl_accINFILE, ctrl_avcOUTNAME, ctrl_avcMSG)
      print*, '... communicated with controller...'
#endif    
      !< Updating the new values
      if (ctrl_aviFAIL .lt. 0) then
        print*, '... error in controller DLL: continuing without controller...'
        this%boolean_abort = .TRUE.
      else
        if (this%swap(1) .eq. 0) print*, '... open controller libdiscon.dll succesfully'
        this%swap       = ctrl_avrSWAP              !< store the new values in the swap array
        this%swap(1)    = 1.0d0
      end if
    end if

    return
end subroutine update_servo

!< subroutine to read input files
subroutine readinginput(this)
    implicit none
    integer :: i, io_error, linejump, UnIn
    character(100) :: str_infile_trim
    class(model_controller), intent(inout)  :: this

    linejump = 3
    
    !< reading contoller_input
    print*, 'reading controller input ...'
    call GetNewUnit (UnIn)
    open(unit = UnIn, file = 'controllerinput.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then 
            do i = 1, linejump
                read(UnIn, *)
            end do
        read(UnIn, *) this%softwaretype
            do i = 1, linejump
                read(UnIn, *)
            end do
        read(UnIn,*) str_infile_trim
            allocate(this%infile, source = trim(adjustl(str_infile_trim)))
            this%length_infile = len(this%infile)
            do i = 1, linejump
                read(UnIn, *)
            end do
        read(UnIn, *) this%init_pitch
            do i = 1, linejump
                read(UnIn, *)
            end do
        read(UnIn, *) this%number_blades
            do i = 1, linejump
                read(UnIn, *) 
            end do
        allocate(this%blade_nbr(this%number_blades))
        allocate(this%revolutej2_ID(this%number_blades))
        allocate(this%rotationl_ID(this%number_blades))
        read(UnIn, *) this%blade_nbr
        read(UnIn, *) this%revolutej2_ID
        read(UnIn, *) this%rotationl_ID        
    else
        this%boolean_abort = .TRUE.
    end if
    close(unit = UnIn)

    return
end subroutine readinginput

!< Subroutine for opening result files to write necessary controller output
  subroutine output_open(this, output_filename, oldsimfilename)
    implicit none
    
    class(model_controller), intent(inout) :: this
    
    integer(8), parameter :: block_size = 2 ** 200
    integer :: i
    
    character(:), allocatable :: output_filename, oldsimfilename
    character(len=1024) :: char_t
    character(len=7) :: str_status = 'replace'
    character(len=6) :: str_position = 'rewind'
    
    if (this%boolean_oldsim) then
      if (output_filename .eq. oldsimfilename) then
        str_status   = 'unknown'
        str_position = 'append'
      end if
    end if
  
    call GetNewUnit(this%UnIn_servofile)
    open(unit = this%UnIn_servofile, file = output_filename // '_servo.dres', action = 'write', status = str_status, POSITION = str_position)
    
    if (.not.allocated(oldsimfilename)) then
      allocate(oldsimfilename, source = 'none')
    endif
    
    if ((this%boolean_oldsim .eqv. .FALSE.) .or. (output_filename .ne. oldsimfilename)) then
      write(this%UnIn_servofile,'(A26)') '!! DeSiO-Controller output'
      write(this%UnIn_servofile,'(A33)') '!! used controller: libdiscon.dll'
      write(this%UnIn_servofile,'(A2)',advance='no') '!!'
      write(this%UnIn_servofile,'(10X,A4,11X)',advance='no') 'time'
      do i = 1,this%number_blades
        write(this%UnIn_servofile,'(10X,A5,I1,9X)',advance='no') 'pitch', i
      end do
      write(this%UnIn_servofile,'(16X,A10,7X,9X,A9,9X,8X,A10,7X,10X,A6,9X,10X,A5,20X,A5,25X,A12,15X,A7,10X,A3,10X)',advance='yes') 'mgenerator', 'omega hub', 'momega hub', 'thrust', 'v_wind', 'm_rotor_aero', 'azimuth', 'P_el'
    end if
    
    return
  end subroutine output_open
  
!< Subroutine for writing result files
  subroutine output_write(this)
    implicit none
    
    class(model_controller), intent(inout) :: this
    
1000 format(1000000f30.15)
1100 format(1000000E27.18E3)
     
    write(this%UnIn_servofile, 1100) this%time, this%pitch_t(1), this%pitch_t(2), this%pitch_t(3), this%mgenerator, this%omega_hub,  this%momega_hub, this%thrust, this%v_wind, this%m_rotor_aero , this%azimuth(1), this%p_el, this%assign_pitch_t(1),this%del_phi(1),this%del_phi(2),this%del_phi(3), this%hub_aux_1(7), this%hub_aux_2(7)

    return
  end subroutine output_write
  
!< Subroutine for closing result files
  subroutine output_close(this)
    implicit none
    class(model_controller), intent(in) :: this
    
    close(this%UnIn_servofile)
    
    return
  end subroutine output_close
  
!< subroutine for displaying discon-controller libdiscon.dll input/output parameter
subroutine ListDisconInfo()
  implicit none
  
   print*, ' swap(1), Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)'
   print*, ' swap(2), Current time (sec) [t in single precision]'
   print*, ' swap(3), Communication interval (sec)'
   print*, ' swap(4), Blade 1 pitch angle (rad) [SrvD input]'
   print*, ' swap(5), Below-rated pitch angle set-point (rad) [SrvD Ptch_SetPnt parameter]'
   print*, ' swap(6), Minimum pitch angle (rad) [SrvD Ptch_Min parameter]'
   print*, ' swap(7), Maximum pitch angle (rad) [SrvD Ptch_Max parameter]'
   print*, ' swap(8), Minimum pitch rate (most negative value allowed) (rad/s) [SrvD PtchRate_Min parameter]'
   print*, ' swap(9), Maximum pitch rate                               (rad/s) [SrvD PtchRate_Max parameter]'
   print*, ' swap(10), 0 = pitch position actuator, 1 = pitch rate actuator (-) [must be 0 for ServoDyn]'
   print*, ' swap(11), Current demanded pitch angle (rad) [I am sending the previous value for blade 1 from the DLL, in the absence of any more information provided in Bladed documentation]'
   print*, ' swap(12), Current demanded pitch rate  (rad/s) [always zero for ServoDyn]'
   print*, ' swap(13), Demanded power (W) [SrvD GenPwr_Dem parameter from input file]'
   print*, ' swap(14), Measured shaft power (W) [SrvD input]'
   print*, ' swap(15), Measured electrical power output (W) [SrvD calculation from previous step; should technically be a state]'
   print*, ' swap(16), Optimal mode gain (Nm/(rad/s)^2) [if torque-speed table look-up not selected in input file, use SrvD Gain_OM parameter, otherwise use 0 (already overwritten in Init routine)]'
   print*, ' swap(17), Minimum generator speed (rad/s) [SrvD GenSpd_MinOM parameter]'
   print*, ' swap(18), Optimal mode maximum speed (rad/s) [SrvD GenSpd_MaxOMp arameter]'
   print*, ' swap(19), Demanded generator speed above rated (rad/s) [SrvD GenSpd_Dem parameter]'
   print*, ' swap(20), Measured generator speed (rad/s) [SrvD input]'
   print*, ' swap(21), Measured rotor speed (rad/s) [SrvD input]'
   print*, ' swap(22), Demanded generator torque above rated (Nm) [SrvD GenTrq_Dem parameter from input file]'
   print*, ' swap(23), Measured generator torque (Nm) [SrvD calculation from previous step; should technically be a state]'
   print*, ' swap(24), Measured yaw error (rad) [SrvD input]'
   !if (dll_data%DLL_NumTrq==0) then  ! Torque-speed table look-up not selected
   print*, ' swap(25), Start of below-rated torque-speed look-up table (Lookup table not in use)'
   !else                 ! Torque-speed table look-up selected
   print*, ' swap(25), Start of below-rated torque-speed look-up table )'
   !endif
   print*, ' swap(26), No. of points in torque-speed look-up table (-) [SrvD DLL_NumTrq parameter]: '
   print*, ' swap(27), Hub wind speed (m/s) [SrvD input]'
   print*, ' swap(28), Pitch control: 0 = collective, 1 = individual (-) [SrvD Ptch_Cntrl parameter]'
   print*, ' swap(29), Yaw control: 0 = yaw rate control, 1 = yaw torque control (-) [must be 0 for ServoDyn] '
   print*, ' swap(30), lade 1 root out-of-plane bending moment (Nm) [SrvD input]'
   print*, ' swap(31), Blade 2 root out-of-plane bending moment (Nm) [SrvD input]'
   print*, ' swap(32), Blade 3 root out-of-plane bending moment (Nm) [SrvD input]'
   print*, ' swap(33), Blade 2 pitch angle (rad) [SrvD input]'
   print*, ' swap(34), Blade 3 pitch angle (rad) [SrvD input]'
   print*, ' swap(37), Nacelle yaw angle from North (rad)'
   print*, ' swap(49), Maximum number of characters in the "MESSAGE" argument (-) [size of ErrMsg argument plus 1 (we add one for the C NULL CHARACTER)]'
   print*, ' swap(50), Number of characters in the "INFILE"  argument (-) [trimmed length of DLL_InFile parameter plus 1 (we add one for the C NULL CHARACTER)]'
   print*, ' swap(51), Number of characters in the "OUTNAME" argument (-) [trimmed length of RootName parameter plus 1 (we add one for the C NULL CHARACTER)]'
   print*, ' swap(53), Tower top fore-aft     acceleration (m/s^2) [SrvD input]'
   print*, ' swap(54), Tower top side-to-side acceleration (m/s^2) [SrvD input]'
   print*, ' swap(60), Rotor azimuth angle (rad) [SrvD input]'
   print*, ' swap(61), Number of blades (-) [SrvD NumBl parameter]'
   print*, ' swap(62), Maximum number of values which can be returned for logging (-)'
   print*, ' swap(63), Record number for start of logging output (-) [set to ]'
   print*, ' swap(64), Maximum number of characters which can be returned in "OUTNAME" (-) [set to (including the C NULL CHARACTER)]'
   print*, ' swap(66), Start of Platform motion -- 1001'
   print*, ' swap(69), Blade 1 root in-plane bending moment (Nm) [SrvD input]'
   print*, ' swap(70), Blade 2 root in-plane bending moment (Nm) [SrvD input]'
   print*, ' swap(71), Blade 3 root in-plane bending moment (Nm) [SrvD input]'
   print*, ' swap(73), Rotating hub My (GL co-ords) (Nm) [SrvD input]'
   print*, ' swap(74), Rotating hub Mz (GL co-ords) (Nm) [SrvD input]'
   print*, ' swap(75), Fixed    hub My (GL co-ords) (Nm) [SrvD input]'
   print*, ' swap(76), Fixed    hub Mz (GL co-ords) (Nm) [SrvD input]'
   print*, ' swap(77), Yaw bearing  My (GL co-ords) (Nm) [SrvD input]'
   print*, ' swap(78), Yaw bearing  Mz (GL co-ords) (Nm) [SrvD input]'
   print*, ' swap(82), Nacelle roll    acceleration (rad/s^2) [SrvD input] -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system'
   print*, ' swap(83), Nacelle nodding acceleration (rad/s^2) [SrvD input] '
   print*, ' swap(84), Nacelle yaw     acceleration (rad/s^2) [SrvD input] -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system'
   print*, ' swap(95), Reserved (SrvD customization: set to SrvD AirDens parameter)'
   print*, ' swap(96), Reserved (SrvD customization: set to SrvD AvgWindSpeed parameter)'
   print*, ' swap(109), Shaft torque (=hub Mx for clockwise rotor) (Nm) [SrvD input]'
   print*, ' swap(110), Thrust - Rotating low-speed shaft force x (GL co-ords) (N) [SrvD input]'
   print*, ' swap(111), Nonrotating low-speed shaft force y (GL co-ords) (N) [SrvD input]'
   print*, ' swap(112), Nonrotating low-speed shaft force z (GL co-ords) (N) [SrvD input]'
   print*, ' swap(117), Controller state [always set to 0]'
   print*, ' swap(129), Maximum extent of the avrSWAP array'
      !   ! Channels with info retrieved from the DLL (from Retrieve_avrSWAP routine)
   print*, ' swap(35), Generator contactor (-) [GenState from previous call to DLL (initialized to 1)]'
   print*, ' swap(36), Shaft brake status (-) [sent to DLL at the next call; anything other than 0 or 1 is an error] '
   print*, ' swap(41), demanded yaw actuator torque [this output is ignored since record 29 is set to 0 by ServoDyn indicating yaw rate control]'
   print*, ' swap(42), demanded individual pitch 1 position'
   print*, ' swap(43), demanded individual pitch 2 position'
   print*, ' swap(44), demanded individual pitch 3 position'
   print*, ' swap(45), Demanded pitch angle (Collective pitch) (rad)'
   print*, ' swap(47), Demanded generator torque (Nm)'
   print*, ' swap(48), Demanded nacelle yaw rate (rad/s)'
   print*, ' swap(55), UNUSED: Pitch override [anything other than 0 is an error in ServoDyn]'
   print*, ' swap(56), UNUSED: Torque override [anything other than 0 is an error in ServoDyn]'
   print*, ' swap(65), Number of variables returned for logging [anything greater than MaxLoggingChannels is an error]'
   print*, ' swap(107), Brake torque demand (used only when avrSWAP(36) is 16)'
   print*, ' swap(120), Airfoil command, blade 1'
   print*, ' swap(121), Airfoil command, blade 2'
   print*, ' swap(122), Airfoil command, blade 3'
   print*, ' swap(63), Number logging channels'

  return
end subroutine ListDisconInfo


end module class_controller
