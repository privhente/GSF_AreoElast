module class_matrix12
use my_math_structure, only: skew, cross, eye
use my_constants_structure, only: i_33
use my_FileIO_functions
use class_beam_element, only: b1matrix, b2matrix, q0matrix, q1matrix, q2matrix, w11matrix, w22matrix 

implicit none
  

!type declarations  
    type :: matrix_12_property
        real(kind = 8)      :: matrix12_matrix(21)
        real(kind = 8)      :: matrix12_m(10) 
        real(kind = 8)      :: alpha_s
        real(kind = 8)      :: alpha_v
        real(kind = 8)      :: alpha_obj
    end type matrix_12_property

    type :: matrix_12_object
        character(len = 9)                  :: matrix12_strain_flag
        integer                             :: node1, node2, nproperty
        real(kind = 8)                      :: fqint_obj(24), fqint(12), fqdyn_AB(24), fv_AB(24)
        real(kind = 8)                      :: kqq(12,12), kqq_obj(24,24), kqv(24,24), kvq(24,24), kvv(24,24), kvv_diss(24,24)   
        real(kind = 8), dimension(24,24)    :: M_tilde, M_tilde_hat
        real(kind = 8)                      :: penergy
        real(kind = 8), dimension(6)        :: s_diss, s_obj_diss, strain1, strain2
        real(kind = 8), dimension(6,6)      :: c_in, c_diss, c_obj_diss
        real(kind = 8), dimension(3, 3)     :: cgg, ckk, cgk
        real(kind = 8), dimension(4,4)      :: c_mass 
        real(kind = 8), dimension(24,24)    :: M_Beam
        real(kind = 8)                      :: deltat
        real(kind = 8)                      :: v_diss(24)
        ! Variablen zu Untersuchung der Lie-Ableitung
        real(kind = 8), dimension(6)        :: lie_u 
        real(kind = 8), dimension(6)        :: num
        real(kind = 8), dimension(6)        :: ua_dot
      
   

    end type matrix_12_object
      
    type :: matrix_12
        type(matrix_12_property), allocatable   :: matrix12_property(:)
        type(matrix_12_object), allocatable     :: matrix12_object(:)
    
    contains
        procedure :: matrix12_lin_stiffness
        procedure :: matrix12_obj_stiffness
        procedure :: matrix12_mass
        procedure :: matrix12_built_mass
        procedure :: matrix12_damping
        procedure :: write_output
        !procedure :: output_open_matrix
    end type matrix_12

    contains

    subroutine matrix12_built_mass(this, i) ! This should be the correct formulation for the Coupling element
    implicit none
    class(matrix_12), intent(inout)                 :: this
    integer, intent(in)                             :: i
    real(kind = 8), dimension(12,12)                :: m 
    
    ! Building the mass and euler inertia tensor from the input  
    this%matrix12_object(i)%c_mass(:,:) = 0.0d0
    this%matrix12_object(i)%c_mass(1,1) = this%matrix12_property(i)%matrix12_m(1)
    this%matrix12_object(i)%c_mass(1,2) = this%matrix12_property(i)%matrix12_m(5)
    this%matrix12_object(i)%c_mass(1,3) = this%matrix12_property(i)%matrix12_m(4)
    this%matrix12_object(i)%c_mass(2,1) = this%matrix12_property(i)%matrix12_m(5)
    this%matrix12_object(i)%c_mass(2,2) = this%matrix12_property(i)%matrix12_m(3)
    this%matrix12_object(i)%c_mass(2,3) = this%matrix12_property(i)%matrix12_m(6)
    this%matrix12_object(i)%c_mass(3,1) = this%matrix12_property(i)%matrix12_m(4)
    this%matrix12_object(i)%c_mass(3,2) = this%matrix12_property(i)%matrix12_m(6)
    this%matrix12_object(i)%c_mass(3,3) = this%matrix12_property(i)%matrix12_m(2)
    
    ! Projection the mass and euler inertia tensor to the DOFs, i_33 = identity matrix 3x3 
    m(1:3,1:3)     = this%matrix12_object(i)%c_mass(1,1)*i_33  
    m(1:3,4:6)     = this%matrix12_object(i)%c_mass(1,2)*i_33
    m(1:3,7:9)     = this%matrix12_object(i)%c_mass(1,3)*i_33
    m(1:3,10:12)   = this%matrix12_object(i)%c_mass(1,4)*i_33
    m(4:6,1:3)     = this%matrix12_object(i)%c_mass(2,1)*i_33
    m(4:6,4:6)     = this%matrix12_object(i)%c_mass(2,2)*i_33
    m(4:6,7:9)     = this%matrix12_object(i)%c_mass(2,3)*i_33
    m(4:6,10:12)   = this%matrix12_object(i)%c_mass(2,4)*i_33
    m(7:9,1:3)     = this%matrix12_object(i)%c_mass(3,1)*i_33
    m(7:9,4:6)     = this%matrix12_object(i)%c_mass(3,2)*i_33
    m(7:9,7:9)     = this%matrix12_object(i)%c_mass(3,3)*i_33
    m(7:9,10:12)   = this%matrix12_object(i)%c_mass(3,4)*i_33
    m(10:12,1:3)   = this%matrix12_object(i)%c_mass(4,1)*i_33
    m(10:12,4:6)   = this%matrix12_object(i)%c_mass(4,2)*i_33
    m(10:12,7:9)   = this%matrix12_object(i)%c_mass(4,3)*i_33
    m(10:12,10:12) = this%matrix12_object(i)%c_mass(4,4)*i_33
    
    ! Assembling M_tilde
    this%matrix12_object(i)%M_tilde(1:12,1:12)   = m
    this%matrix12_object(i)%M_tilde(1:12,13:24)  = m
    this%matrix12_object(i)%M_tilde(13:24,1:12)  = m
    this%matrix12_object(i)%M_tilde(13:24,13:24) = m
    this%matrix12_object(i)%M_tilde(:,:)         = 0.25d0*this%matrix12_object(i)%M_tilde
    
    this%matrix12_object(i)%M_tilde_hat(10:12,10:12) = eye(3)
    this%matrix12_object(i)%M_tilde_hat(22:24,22:24) = eye(3)
    this%matrix12_object(i)%M_tilde_hat(10:12,22:24) = eye(3)
    this%matrix12_object(i)%M_tilde_hat(22:24,10:12) = eye(3)
    
    
    this%matrix12_object(i)%M_Beam =  beammassmatrix_CE(this) ! using copied beam element function to compare to beam mass matrix

        
    end subroutine matrix12_built_mass

    
!----------------------- Functions to bulit beam mass matrix ------------------------------------------------------------------------
    
    !function beammassmatrix_CE(this,q0a,q0b) result(M_Beam) 
    function beammassmatrix_CE(this) result(M_Beam) 
    use class_beam_element
    implicit none

    class(matrix_12), intent(inout)  :: this
    
    real(kind = 8), dimension(24, 24) :: massmatrix
    real(kind = 8), dimension(24, 24) :: mhat

    integer :: i
    real(kind = 8), dimension(3)     :: deltaphi_0
    real(kind = 8)                   :: length_0
    real(kind = 8)                   :: jacobian
    real(kind = 8), dimension(2)     :: n1, n2, weigth_m
    real(kind = 8), dimension(24)    :: q_0
    real(kind = 8), dimension(12)    :: q0a, q0b
    real(kind = 8), dimension(12,24) :: q0
    real(kind = 8), dimension(12,12) :: w00
    real(kind = 8), dimension(24,24) :: m, M_Beam

    
    !deltaphi_0 = q0a(1:3)-q0b(1:3)   
    !length_0 = dsqrt(dot_product(deltaphi_0, deltaphi_0))
    length_0 = 1.0d0
    
    n1(1) = 0.788675134594813d0
    n1(2) = 0.211324865405187d0
    n2(1) = n1(2)
    n2(2) = n1(1)
    jacobian = 0.5d0*length_0
    weigth_m(:) = 1.0d0
    
    w00(:, :) = 0.0d0

    w00(1:3, 1:3) = this%matrix12_property(1)%matrix12_m(1)*i_33
    w00(1:3, 4:6) = this%matrix12_property(1)%matrix12_m(5)*i_33
    w00(1:3, 7:9) = this%matrix12_property(1)%matrix12_m(4)*i_33

    w00(4:6, 1:3) = this%matrix12_property(1)%matrix12_m(5)*i_33
    w00(4:6, 4:6) = this%matrix12_property(1)%matrix12_m(3)*i_33
    w00(4:6, 7:9) = this%matrix12_property(1)%matrix12_m(6)*i_33
    
    w00(7:9, 1:3) = this%matrix12_property(1)%matrix12_m(4)*i_33
    w00(7:9, 4:6) = this%matrix12_property(1)%matrix12_m(6)*i_33
    w00(7:9, 7:9) = this%matrix12_property(1)%matrix12_m(2)*i_33
  
    
    massmatrix(:, :) = 0.0d0

    m(:, :)   = 0.0d0
    mhat(:,:) = 0.0d0
    
    ! 2 point Gauss integration 
    do i = 1, 2 

       q0 = q0_matrix(n1(i), n2(i))
       
       massmatrix = massmatrix+weigth_m(i)*matmul(transpose(q0), matmul(w00, q0))
       
    end do

    massmatrix =  massmatrix*jacobian 
    M_Beam(:,:) = massmatrix
 
    return
    
    end function beammassmatrix_CE
    
    function q0_matrix(n1, n2)
    implicit none
    real(kind = 8), dimension(12, 24) :: q0_matrix
    real(kind = 8), intent(in) :: n1, n2

    q0_matrix(:, :) = 0.0d0

    q0_matrix( 1: 3,  1: 3) = n1*i_33
    q0_matrix( 4: 6,  4: 6) = n1*i_33
    q0_matrix( 7: 9,  7: 9) = n1*i_33
    q0_matrix(10:12, 10:12) = n1*i_33

    q0_matrix( 1: 3, 13:15) = n2*i_33
    q0_matrix( 4: 6, 16:18) = n2*i_33
    q0_matrix( 7: 9, 19:21) = n2*i_33
    q0_matrix(10:12, 22:24) = n2*i_33

    return
    end function q0_matrix 
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
    
subroutine matrix12_lin_stiffness(this, qtn1, qtn, q0, number_matrix12, i)
    implicit none
        
    class(matrix_12), intent(inout)  :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    integer, intent(in)              :: number_matrix12, i
        
    real(kind = 8), dimension(12,12) :: kg1, kg2, km
    real(kind = 8), dimension(6,12)  :: d_matrix1
    real(kind = 8), dimension(6,12)  :: d_matrix2
    real(kind = 8), dimension(6)     :: smatrix, b
    real(kind = 8), dimension(6,6)   :: m_in
    real(kind = 8), dimension(12)    :: qt05
     
        ! Building Matrix C_in
        m_in(:,:) = 0.0d0
        m_in(1,2) = this%matrix12_property(i)%matrix12_matrix(15)
        m_in(1,3) = this%matrix12_property(i)%matrix12_matrix(14)
        m_in(1,4) = this%matrix12_property(i)%matrix12_matrix(13)
        m_in(1,5) = this%matrix12_property(i)%matrix12_matrix(12)
        m_in(1,6) = this%matrix12_property(i)%matrix12_matrix(11)
        m_in(2,3) = this%matrix12_property(i)%matrix12_matrix(16)
        m_in(2,4) = this%matrix12_property(i)%matrix12_matrix(21)
        m_in(2,5) = this%matrix12_property(i)%matrix12_matrix(20)
        m_in(2,6) = this%matrix12_property(i)%matrix12_matrix(10)
        m_in(3,4) = this%matrix12_property(i)%matrix12_matrix(17)
        m_in(3,5) = this%matrix12_property(i)%matrix12_matrix(19)
        m_in(3,6) = this%matrix12_property(i)%matrix12_matrix(9)
        m_in(4,5) = this%matrix12_property(i)%matrix12_matrix(18)
        m_in(4,6) = this%matrix12_property(i)%matrix12_matrix(8)       
        m_in(5,6) = this%matrix12_property(i)%matrix12_matrix(7)
          
        m_in = m_in+transpose(m_in)
        m_in(1,1) = this%matrix12_property(i)%matrix12_matrix(1)
        m_in(2,2) = this%matrix12_property(i)%matrix12_matrix(2)
        m_in(3,3) = this%matrix12_property(i)%matrix12_matrix(3)
        m_in(4,4) = this%matrix12_property(i)%matrix12_matrix(4)
        m_in(5,5) = this%matrix12_property(i)%matrix12_matrix(5)
        m_in(6,6) = this%matrix12_property(i)%matrix12_matrix(6)
        
        ! Displacement at t = 0.5 
        qt05(:) = 0.5d0*(qtn1(:) + qtn(:)) 
        
        ! Defining the translation and rotation in global coordinates
        b(1:3) = qt05(1:3) - q0(1:3)
        b(4:6) = -0.5d0*(cross(qt05(4:6),q0(4:6))+cross(qt05(7:9),q0(7:9))+cross(qt05(10:12),q0(10:12)))  
          
        !! Calculation S (Vector of Forces and Moments)
        smatrix = matmul(m_in, b)
        
        ! Calculation D1
        d_matrix1(:,:) = 0.0d0
        d_matrix1(1:3,1:3)   = eye(3)
        d_matrix1(4:6,4:6)   = 0.5d0*skew(qt05(4:6))
        d_matrix1(4:6,7:9)   = 0.5d0*skew(qt05(7:9))
        d_matrix1(4:6,10:12) = 0.5d0*skew(qt05(10:12))
        
        ! Calculation D2
        d_matrix2(:,:)         =  0.0d0
        d_matrix2(4:6,4:6)     = -0.5d0*skew(qt05(4:6)-q0(4:6))
        d_matrix2(4:6,7:9)     = -0.5d0*skew(qt05(7:9)-q0(7:9))
        d_matrix2(4:6,10:12)   = -0.5d0*skew(qt05(10:12)-q0(10:12))
 
        ! Calculation kg1
        kg1(:,:)               =  0.0d0
        kg1(4:6,4:6)           = -0.5d0*skew(smatrix(4:6))
        kg1(7:9,7:9)           = kg1(4:6,4:6)
        kg1(10:12,10:12)       = kg1(4:6,4:6)
          
        ! Calculation kg2
        kg2 = matmul(matmul(transpose(d_matrix1),m_in),d_matrix2)
          
        ! Calculation km (material stiffness) 
        km = matmul(matmul(transpose(d_matrix1),m_in),d_matrix1)
        
        ! Added Stiffness
        this%matrix12_object(i)%kqq = 0.5d0*(kg1+kg2+km)     
        
        ! Calculation fqint 
        this%matrix12_object(i)%fqint(:) = 0.0d0
        this%matrix12_object(i)%fqint(:) = matmul(transpose(d_matrix1),smatrix)        
    
        ! Calculation Energy
        this%matrix12_object(i)%penergy  = 0.5d0*dot_product(qtn1-q0,matmul(this%matrix12_object(i)%kqq, qtn1-q0)) ! E = 0,5*k*x^2          
          
return
end subroutine matrix12_lin_stiffness

!------------------------------------------------------------------------------------------------------------------------------
subroutine matrix12_obj_stiffness(this, qtn1a_opt, qtna, q0a, qtn1b_opt, qtnb, q0b, vtn1a_opt, vtna_opt, vtn1b_opt, vtnb_opt, number_matrix12, i)
    
implicit none 
        
    class(matrix_12), intent(inout)                      :: this
    real(kind = 8), dimension(12), intent(in)            :: qtna, q0a, qtnb, q0b
    real(kind = 8), dimension(12), optional,  intent(in) :: qtn1a_opt, qtn1b_opt
    integer, intent(in)                                  :: number_matrix12, i
    real(kind = 8), dimension(12), optional, intent(in)  :: vtn1a_opt, vtna_opt, vtn1b_opt, vtnb_opt               ! Velocities
    integer, dimension(24)                               :: indices_q_total
        
    real(kind = 8), dimension(24,24) :: kg1, kg2, km
    real(kind = 8), dimension(12)    :: qt05a, qt05b, qtn1a, qtn1b
    real(kind = 8), dimension(3)     :: d21a, d22a, d23a, phi2a         ! displacement node A, time t+1 
    real(kind = 8), dimension(3)     :: d21b, d22b, d23b, phi2b         ! displacement node B, time t+1
    real(kind = 8), dimension(3)     :: d11a, d12a, d13a, phi1a         ! displacement node A, time t 
    real(kind = 8), dimension(3)     :: d11b, d12b, d13b, phi1b         ! displacement node B, time t
    real(kind = 8), dimension(3)     :: d01a, d02a, d03a, phi0a         ! displacement node A, time 0 
    real(kind = 8), dimension(3)     :: d01b, d02b, d03b, phi0b         ! displacement node B, time 0
    real(kind = 8), dimension(3)     :: d1a05, d2a05, d3a05, phia05     ! displacement node A, time 0.5 
    real(kind = 8), dimension(3)     :: d1b05, d2b05, d3b05, phib05     ! displacement node B, time 0.5
    real(kind = 8), dimension(6,24)  :: b1_matrix05, b1_matrix2, b2_matrix
    real(kind = 8), dimension(6)     :: nmatrix_05, nmatrix1, nmatrix2, u1, u2
    real(kind = 8), dimension(6,6)   :: m_in
    real(kind = 8)                   :: eta
    real(kind = 8)                   :: length0, jacobian, weight
    real(kind = 8), dimension(12)    :: vtn1a, vtna, vtn1b, vtnb               ! velocities

    
    if (present(qtn1a_opt)) then
       qtn1a = qtn1a_opt
    else
       qtn1a = qtna
    end if
    
    if (present(qtn1b_opt)) then
       qtn1b = qtn1b_opt
    else
       qtn1b = qtnb
    end if

    if (present(vtna_opt) .and. present(vtn1a_opt)) then
       vtna = vtna_opt
       vtn1a = vtn1a_opt
    else
       vtna(:) = 0.0d0
       vtn1a(:) = 0.0d0
    end if
    
    if (present(vtnb_opt) .and. present(vtn1b_opt)) then
       vtnb = vtnb_opt
       vtn1b = vtn1b_opt
    else
       vtna(:) = 0.0d0
       vtn1b(:) = 0.0d0
    end if
    ! -----------------------------------------  
    ! Defining the directors at different times
    ! -----------------------------------------  

        ! t0, node A 
        phi0a  = q0a(1:3)
        d01a   = q0a(4:6)
        d02a   = q0a(7:9)
        d03a   = q0a(10:12)
    
        ! tn, node A
        phi1a  = qtna(1:3)
        d11a   = qtna(4:6)
        d12a   = qtna(7:9)
        d13a   = qtna(10:12)    
    
        ! tn+1, node A
        phi2a  = qtn1a(1:3)
        d21a   = qtn1a(4:6)
        d22a   = qtn1a(7:9)
        d23a   = qtn1a(10:12) 
    
        ! t0, node B
        phi0b  = q0b(1:3)
        d01b   = q0b(4:6)
        d02b   = q0b(7:9)
        d03b   = q0b(10:12)
    
        ! tn, node B
        phi1b  = qtnb(1:3)
        d11b   = qtnb(4:6)
        d12b   = qtnb(7:9)
        d13b   = qtnb(10:12)    
    
        ! tn+1, node B
        phi2b  = qtn1b(1:3)
        d21b   = qtn1b(4:6)
        d22b   = qtn1b(7:9)
        d23b   = qtn1b(10:12) 
        
        ! Displacement A  at tn+1/2
        qt05a(:) = 0.5d0*(qtn1a(:) + qtna(:)) 

        ! Defining A directors at t = 0.5
        phia05  = qt05a(1:3)
        d1a05   = qt05a(4:6)
        d2a05   = qt05a(7:9)
        d3a05   = qt05a(10:12)
        
        ! Displacement B at t = 0.5
        qt05b(:) = 0.5d0*(qtn1b(:) + qtnb(:)) 

        ! Directors B at t = 0.5
        phib05  = qt05b(1:3)
        d1b05   = qt05b(4:6)
        d2b05   = qt05b(7:9)
        d3b05   = qt05b(10:12)

        
    ! -----------------------------------------  
    ! Parametrization for the beam formulation
    ! -----------------------------------------  

        length0  = dsqrt(dot_product(q0b-q0a,q0b-q0a))
        jacobian = 0.5d0*length0
        weight   = 2.0d0
        eta      = length0    
        
        
        !General couling element formulation
        !length0  = 1.0d0
        !jacobian = 1.0d0
        !weight   = 1.0d0
        !eta      = 1.0d0   
       
        
    
    ! -----------------------------------------  
    ! Building Matrix C_in --> 11 Parameter version   !TODO: built in flag
    ! -----------------------------------------  
       m_in(:,:) = 0.0d0
       m_in(1,6) = -this%matrix12_property(i)%matrix12_matrix(9)
       m_in(2,6) = this%matrix12_property(i)%matrix12_matrix(10)
       m_in(3,4) = this%matrix12_property(i)%matrix12_matrix(7)
       m_in(3,5) = -this%matrix12_property(i)%matrix12_matrix(8)
       m_in(4,5) = -this%matrix12_property(i)%matrix12_matrix(11)
       
       m_in(6,1) = -this%matrix12_property(i)%matrix12_matrix(9)
       m_in(6,2) = this%matrix12_property(i)%matrix12_matrix(10)
       m_in(4,3) = this%matrix12_property(i)%matrix12_matrix(7)
       m_in(5,3) = -this%matrix12_property(i)%matrix12_matrix(8)
       m_in(5,4) = -this%matrix12_property(i)%matrix12_matrix(11)
       
       m_in(1,1) = this%matrix12_property(i)%matrix12_matrix(2)
       m_in(2,2) = this%matrix12_property(i)%matrix12_matrix(3)
       m_in(3,3) = this%matrix12_property(i)%matrix12_matrix(1)
       m_in(4,4) = this%matrix12_property(i)%matrix12_matrix(4)
       m_in(5,5) = this%matrix12_property(i)%matrix12_matrix(5)
       m_in(6,6) = this%matrix12_property(i)%matrix12_matrix(6)
       
       
    ! -----------------------------------------  
    ! Building Matrix C_in --> Voigt Notation
    ! -----------------------------------------  
       !m_in(:,:) = 0.0d0
       !
       !! top right
       !m_in(1,2) = -this%matrix12_property(i)%matrix12_matrix(15)
       !m_in(1,3) = -this%matrix12_property(i)%matrix12_matrix(14)
       !m_in(1,4) = -this%matrix12_property(i)%matrix12_matrix(13)
       !m_in(1,5) = -this%matrix12_property(i)%matrix12_matrix(12)
       !m_in(1,6) = -this%matrix12_property(i)%matrix12_matrix(11)
       !            
       !m_in(2,3) = -this%matrix12_property(i)%matrix12_matrix(16)
       !m_in(2,4) = -this%matrix12_property(i)%matrix12_matrix(21)
       !m_in(2,5) = -this%matrix12_property(i)%matrix12_matrix(20)
       !m_in(2,6) = -this%matrix12_property(i)%matrix12_matrix(10)
       !            
       !m_in(3,4) = -this%matrix12_property(i)%matrix12_matrix(17)
       !m_in(3,5) = -this%matrix12_property(i)%matrix12_matrix(19)
       !m_in(3,6) = -this%matrix12_property(i)%matrix12_matrix(9)
       !            
       !m_in(4,5) = -this%matrix12_property(i)%matrix12_matrix(18)
       !m_in(4,6) = -this%matrix12_property(i)%matrix12_matrix(8)
       !            
       !m_in(5,6) = -this%matrix12_property(i)%matrix12_matrix(7)
       !
       !! bottom left
       !m_in(2,1) = -this%matrix12_property(i)%matrix12_matrix(15)
       !m_in(3,1) = -this%matrix12_property(i)%matrix12_matrix(14)
       !m_in(4,1) = -this%matrix12_property(i)%matrix12_matrix(13)
       !m_in(5,1) = -this%matrix12_property(i)%matrix12_matrix(12)
       !m_in(6,1) = -this%matrix12_property(i)%matrix12_matrix(11)
       !            
       !m_in(3,2) = -this%matrix12_property(i)%matrix12_matrix(16)
       !m_in(4,2) = -this%matrix12_property(i)%matrix12_matrix(21)
       !m_in(5,2) = -this%matrix12_property(i)%matrix12_matrix(20)
       !m_in(6,2) = -this%matrix12_property(i)%matrix12_matrix(10)
       !            
       !m_in(4,3) = -this%matrix12_property(i)%matrix12_matrix(17)
       !m_in(5,3) = -this%matrix12_property(i)%matrix12_matrix(19)
       !m_in(6,3) = -this%matrix12_property(i)%matrix12_matrix(9)
       !            
       !m_in(5,4) = -this%matrix12_property(i)%matrix12_matrix(18)
       !m_in(6,4) = -this%matrix12_property(i)%matrix12_matrix(8)
       !           
       !m_in(6,5) = -this%matrix12_property(i)%matrix12_matrix(7)
       !
       !! diagonal
       !m_in(1,1) = this%matrix12_property(i)%matrix12_matrix(1)
       !m_in(2,2) = this%matrix12_property(i)%matrix12_matrix(2)
       !m_in(3,3) = this%matrix12_property(i)%matrix12_matrix(3)
       !m_in(4,4) = this%matrix12_property(i)%matrix12_matrix(4)
       !m_in(5,5) = this%matrix12_property(i)%matrix12_matrix(5)
       !m_in(6,6) = this%matrix12_property(i)%matrix12_matrix(6)
       
    !< store in type matrix to have acces in the whole class. Otherwise rename m_in to this%---%c_in
       this%matrix12_object(i)%c_in = m_in  
       
    ! -----------------------------------------  
    ! Defining the objective nonlinear deformation measure u = [v w]^T / [Gamma, Kappa]^T
    ! -----------------------------------------  
        ! deformation/strain at tn
        u1(1) = (1.0d0/(2.0d0 * eta)) * (dot_product((phi1b-phi1a),(d11a+d11b)) - dot_product((phi0b-phi0a),(d01a+d01b)))
        u1(2) = (1.0d0/(2.0d0 * eta)) * (dot_product((phi1b-phi1a),(d12a+d12b)) - dot_product((phi0b-phi0a),(d02a+d02b)))
        u1(3) = (1.0d0/(2.0d0 * eta)) * (dot_product((phi1b-phi1a),(d13a+d13b)) - dot_product((phi0b-phi0a),(d03a+d03b)))
        u1(4) = (1.0d0/(2.0d0 * eta)) * (dot_product(d12a, d13b) - dot_product(d13a, d12b) -  (dot_product(d02a, d03b) - dot_product(d03a, d02b)))
        u1(5) = (1.0d0/(2.0d0 * eta)) * (dot_product(d13a, d11b) - dot_product(d11a, d13b) -  (dot_product(d03a, d01b) - dot_product(d01a, d03b)))
        u1(6) = (1.0d0/(2.0d0 * eta)) * (dot_product(d11a, d12b) - dot_product(d12a, d11b) -  (dot_product(d01a, d02b) - dot_product(d02a, d01b)))
        
        ! deformation/strain at tn+1
        u2(1) = (1.0d0/(2.0d0 * eta)) * (dot_product((phi2b-phi2a),(d21a+d21b)) - dot_product((phi0b-phi0a),(d01a+d01b)))
        u2(2) = (1.0d0/(2.0d0 * eta)) * (dot_product((phi2b-phi2a),(d22a+d22b)) - dot_product((phi0b-phi0a),(d02a+d02b)))
        u2(3) = (1.0d0/(2.0d0 * eta)) * (dot_product((phi2b-phi2a),(d23a+d23b)) - dot_product((phi0b-phi0a),(d03a+d03b)))
        u2(4) = (1.0d0/(2.0d0 * eta)) * (dot_product(d22a, d23b) - dot_product(d23a, d22b) -  (dot_product(d02a, d03b) - dot_product(d03a, d02b)))
        u2(5) = (1.0d0/(2.0d0 * eta)) * (dot_product(d23a, d21b) - dot_product(d21a, d23b) -  (dot_product(d03a, d01b) - dot_product(d01a, d03b)))
        u2(6) = (1.0d0/(2.0d0 * eta)) * (dot_product(d21a, d22b) - dot_product(d22a, d21b) -  (dot_product(d01a, d02b) - dot_product(d02a, d01b)))
        
        this%matrix12_object(i)%strain1 = u1
        this%matrix12_object(i)%strain2 = u2
    ! -----------------------------------------  
    ! Calculation N (Vector of Forces and Moments)
    ! -----------------------------------------  
        ! Calculation of the dissipation
        call this%matrix12_damping(qtn1a, qtna, qtn1b, qtnb, vtn1a, vtn1b, vtna, vtnb, i)
                           
        ! Calculation N 
        nmatrix1 = 0.0d0
        nmatrix1 = matmul(m_in, u1)            ! Forces and Moments at tn
        
        nmatrix2 = 0.0d0
        nmatrix2 = matmul(m_in, u2)            ! Forces an Moments at tn+1
        
        nmatrix_05 = 0.0d0
        nmatrix_05 = 0.5d0*(nmatrix1+nmatrix2) ! Average forces and moments between tn and tn+1
        
        nmatrix_05 = nmatrix_05 + this%matrix12_object(i)%s_diss + this%matrix12_object(i)%s_obj_diss
        
    ! -----------------------------------------  
    ! Calculation B1 at tn+1/2
    ! -----------------------------------------  
        b1_matrix05(:,:)      = 0.0d0
        
        b1_matrix05(1,1:3)    = -d1a05-d1b05
        b1_matrix05(1,4:6)    =  phib05-phia05
        b1_matrix05(1,13:15)  =  d1a05+d1b05
        b1_matrix05(1,16:18)  =  phib05-phia05
        
        b1_matrix05(2,1:3)    = -d2a05-d2b05
        b1_matrix05(2,7:9)    =  phib05-phia05
        b1_matrix05(2,13:15)  =  d2a05+d2b05
        b1_matrix05(2,19:21)  =  phib05-phia05
        
        b1_matrix05(3,1:3)    = -d3a05-d3b05
        b1_matrix05(3,10:12)  =  phib05-phia05
        b1_matrix05(3,13:15)  =  d3a05+d3b05
        b1_matrix05(3,22:24)  =  phib05-phia05
     
        b1_matrix05(4,7:9)    =  d3b05
        b1_matrix05(4,10:12)  = -d2b05
        b1_matrix05(4,19:21)  = -d3a05
        b1_matrix05(4,22:24)  =  d2a05
  
        b1_matrix05(5,4:6)    = -d3b05
        b1_matrix05(5,10:12)  =  d1b05  
        b1_matrix05(5,16:18)  =  d3a05
        b1_matrix05(5,22:24)  = -d1a05
    
        b1_matrix05(6,4:6)    =  d2b05
        b1_matrix05(6,7:9)    = -d1b05 
        b1_matrix05(6,16:18)  = -d2a05
        b1_matrix05(6,19:21)  =  d1a05
        
        b1_matrix05 = 1.0d0/(2.0d0*eta) * b1_matrix05
        
    ! -----------------------------------------  
    ! Calculation B1 at tn+1
    ! -----------------------------------------  
        b1_matrix2(:,:)      = 0.0d0
        
        b1_matrix2(1,1:3)    = -d21a-d21b
        b1_matrix2(1,4:6)    =  phi2b-phi2a
        b1_matrix2(1,13:15)  =  d21a+d21b
        b1_matrix2(1,16:18)  =  phi2b-phi2a
        
        b1_matrix2(2,1:3)    = -d22a-d22b
        b1_matrix2(2,7:9)    =  phi2b-phi2a
        b1_matrix2(2,13:15)  =  d22a+d22b
        b1_matrix2(2,19:21)  =  phi2b-phi2a
        
        b1_matrix2(3,1:3)    = -d23a-d23b
        b1_matrix2(3,10:12)  =  phi2b-phi2a
        b1_matrix2(3,13:15)  =  d23a+d23b
        b1_matrix2(3,22:24)  =  phi2b-phi2a
     
        b1_matrix2(4,7:9)    =  d23b
        b1_matrix2(4,10:12)  = -d22b
        b1_matrix2(4,19:21)  = -d23a
        b1_matrix2(4,22:24)  =  d22a
  
        b1_matrix2(5,4:6)    = -d23b
        b1_matrix2(5,10:12)  =  d21b  
        b1_matrix2(5,16:18)  =  d23a
        b1_matrix2(5,22:24)  = -d21a
    
        b1_matrix2(6,4:6)    =  d22b
        b1_matrix2(6,7:9)    = -d21b 
        b1_matrix2(6,16:18)  = -d22a
        b1_matrix2(6,19:21)  =  d21a
        
        b1_matrix2 = 1.0d0/(2.0d0*eta) * b1_matrix2
        
    ! -----------------------------------------  
    ! Calculation B2 at tn+1/2
    ! -----------------------------------------  
        b2_matrix(:,:)      = 0.0d0
        
        b2_matrix(1,1:3)    = -(d1a05-d01a) - (d1b05- d01b)
        b2_matrix(1,4:6)    =  (phib05-phi0b)-(phia05-phi0a)
        b2_matrix(1,13:15)  =  (d1a05-d1a05) +(d1b05- d01b)
        b2_matrix(1,16:18)  =  (phib05-phi0b)-(phia05-phi0a)
        
        b2_matrix(2,1:3)    = -(d2a05-d02a) - (d2b05 -d02b)
        b2_matrix(2,7:9)    =  (phib05-phi0b)-(phia05-phi0a)
        b2_matrix(2,13:15)  =  (d2a05-d02a)  +(d2b05 -d02b)
        b2_matrix(2,19:21)  =  (phib05-phi0b)-(phia05-phi0a)
                            
        b2_matrix(3,1:3)    = -(d3a05-d03a)  -(d3b05 -d03b)
        b2_matrix(3,10:12)  =  (phib05-phi0b)-(phia05-phi0a)
        b2_matrix(3,13:15)  =  (d3a05-d03a)  +(d3b05 -d03b)
        b2_matrix(3,22:24)  =  (phib05-phi0b)-(phia05-phi0a)
     
        b2_matrix(4,7:9)    =  (d3b05-d03b)
        b2_matrix(4,10:12)  = -(d2b05-d02b) 
        b2_matrix(4,19:21)  = -(d3a05-d03a)
        b2_matrix(4,22:24)  =  (d2a05-d02a)
  
        b2_matrix(5,4:6)    = -(d3b05-d03b)
        b2_matrix(5,10:12)  =  (d1b05-d01b) 
        b2_matrix(5,16:18)  =  (d3a05-d03a)
        b2_matrix(5,22:24)  = -(d1a05-d01a)
    
        b2_matrix(6,4:6)    =  (d2b05-d02b)
        b2_matrix(6,7:9)    = -(d1b05-d01b)
        b2_matrix(6,16:18)  = -(d2a05-d02a)
        b2_matrix(6,19:21)  =  (d1a05-d01a)
        
        b2_matrix = 1.0d0/(2.0d0*eta) * b2_matrix
    
        ! Calculation kg1
        kg1(:,:) = 0.0d0
        
        kg1(1:3,4:6)      = -nmatrix_05(1)*eye(3) 
        kg1(1:3,7:9)      = -nmatrix_05(2)*eye(3) 
        kg1(1:3,10:12)    = -nmatrix_05(3)*eye(3) 
        kg1(1:3,16:18)    = -nmatrix_05(1)*eye(3) 
        kg1(1:3,19:21)    = -nmatrix_05(2)*eye(3) 
        kg1(1:3,22:24)    = -nmatrix_05(3)*eye(3)
        
        kg1(4:6,1:3)      = -nmatrix_05(1)*eye(3) 
        kg1(4:6,13:15)    = nmatrix_05(1)*eye(3) 
        kg1(4:6,19:21)    = nmatrix_05(6)*eye(3) 
        kg1(4:6,22:24)    = -nmatrix_05(5)*eye(3) 
        
        kg1(7:9,1:3)      = -nmatrix_05(2)*eye(3)
        kg1(7:9,13:15)    = nmatrix_05(2)*eye(3)
        kg1(7:9,16:18)    = -nmatrix_05(6)*eye(3)
        kg1(7:9,22:24)    = nmatrix_05(4)*eye(3)
        
        kg1(10:12,1:3)    = -nmatrix_05(3)*eye(3)
        kg1(10:12,13:15)  = nmatrix_05(3)*eye(3)
        kg1(10:12,16:18)  = nmatrix_05(5)*eye(3)
        kg1(10:12,19:21)  = -nmatrix_05(4)*eye(3)
        
        kg1(13:15,4:6)    = nmatrix_05(1)*eye(3)
        kg1(13:15,7:9)    = nmatrix_05(2)*eye(3)
        kg1(13:15,10:12)  = nmatrix_05(3)*eye(3)
        kg1(13:15,16:18)  = nmatrix_05(1)*eye(3)
        kg1(13:15,19:21)  = nmatrix_05(2)*eye(3)
        kg1(13:15,22:24)  = nmatrix_05(3)*eye(3)
        
        kg1(16:18,1:3)    = -nmatrix_05(1)*eye(3)
        kg1(16:18,7:9)    = -nmatrix_05(6)*eye(3)
        kg1(16:18,10:12)  = nmatrix_05(5)*eye(3)
        kg1(16:18,13:15)  = nmatrix_05(1)*eye(3)
        
        kg1(19:21,1:3)    = -nmatrix_05(2)*eye(3)
        kg1(19:21,4:6)    = nmatrix_05(6)*eye(3)
        kg1(19:21,10:12)  = -nmatrix_05(4)*eye(3)
        kg1(19:21,13:15)  = nmatrix_05(2)*eye(3)
        
        kg1(22:24,1:3)    = -nmatrix_05(3)*eye(3)
        kg1(22:24,4:6)    = -nmatrix_05(5)*eye(3)
        kg1(22:24,7:9)    = nmatrix_05(4)*eye(3)
        kg1(22:24,13:15)  = nmatrix_05(3)*eye(3)

    !------------------------
    ! Calculation of the internal force and the tangential matrix
    !------------------------
        ! Geometric stiffness
        kg1 = (1.0d0/(2.0d0*eta))*kg1
        kg2(:,:) = 0.0d0
        !kg2 = matmul(transpose(b1_matrix05), matmul(m_in, b2_matrix))
        
        ! Material stiffness 
        km(:,:) = 0.0d0
        km(:,:) = matmul(transpose(b1_matrix05), matmul(m_in+2*this%matrix12_object(i)%c_diss+2*this%matrix12_object(i)%c_obj_diss, b1_matrix2))
             
        ! Added Stiffness
        this%matrix12_object(i)%kqq_obj(:,:) = 0.0d0
        this%matrix12_object(i)%kqq_obj(:,:) = 0.5d0*(kg1+kg2+km)  
        this%matrix12_object(i)%kqq_obj = weight*jacobian*this%matrix12_object(i)%kqq_obj

        ! Calculation fqint 
        this%matrix12_object(i)%fqint_obj(:)  = 0.0d0 
        this%matrix12_object(i)%fqint_obj(:)  = weight*jacobian* matmul(transpose(b1_matrix05), nmatrix_05) 
        
    !------------------------
    ! Updating the Energy
    !------------------------
        this%matrix12_object(i)%penergy  = 0.5d0*weight*jacobian*dot_product(u2, nmatrix2)
        !this%matrix12_object(i)%penergy  = 0.5d0*dot_product(qtn1a-q0a,matmul(this%matrix12_object(i)%kqq, qtn1a-q0a))+0.5d0*dot_product(qtn1b-q0b,matmul(this%matrix12_object(i)%kqq, qtn1b-q0b))
        
        
return
end subroutine matrix12_obj_stiffness
!------------------------------------------------------------------------------------------------------------------------------
subroutine matrix12_mass(this, deltat, qtn1a, qtna, q0a, qtn1b, qtnb, q0b, vtn1a_opt, vtna_opt, vtn1b_opt, vtnb_opt, number_matrix12, i)
  implicit none
        
    class(matrix_12), intent(inout)             :: this
    real(kind = 8), dimension(12), intent(in)   :: qtn1a, qtna, q0a, qtn1b, qtnb, q0b     ! Positions
    real(kind = 8), dimension(12)               :: vtn1a, vtna, vtn1b, vtnb               ! Velocities
    real(kind = 8), dimension(12), optional, intent(in)   :: vtn1a_opt, vtna_opt, vtn1b_opt, vtnb_opt               ! Velocities
    real(kind = 8), intent(in)                  :: deltat
    integer, intent(in)                         :: number_matrix12, i 
    real(kind = 8), dimension(24)               :: qtn1AB, qtnAB, vtn1AB, vtnAB, aAB  , vq, vav    ! Pos., Veloc. and Acc. for point C
    real(kind = 8), dimension(12, 12)           :: m
    real(kind = 8), dimension(24,24)            :: M_Beam_hat

    if (present(vtna_opt) .and. present(vtn1a_opt)) then
       vtna = vtna_opt
       vtn1a = vtn1a_opt
    else
       vtna(:) = 0.0d0
       vtn1a(:) = 0.0d0
    end if
    
    if (present(vtnb_opt) .and. present(vtn1b_opt)) then
       vtnb = vtnb_opt
       vtn1b = vtn1b_opt
    else
       vtna(:) = 0.0d0
       vtn1b(:) = 0.0d0
    end if
    
! Calculate position, velocity and acceleration as a single vector for both nodes
      
    qtn1AB(1:12)  = qtn1a(:)
    qtn1AB(13:24) = qtn1b(:)
    qtnAB(1:12)   = qtna(:)
    qtnAB(13:24)  = qtnb(:)
    
    vtn1AB(1:12)  = vtn1a(:)
    vtn1AB(13:24) = vtn1b(:)
    vtnAB(1:12)   = vtna(:)
    vtnAB(13:24)  = vtnb(:)
    
    aAB(:)  = 0.0d0
    aAB(:)  = (vtn1AB(:)-vtnAB(:))/deltat 
    
    ! Beam
    M_Beam_hat(:,:) = this%matrix12_object(i)%M_Beam(:,:)    
    M_Beam_hat(10:12,10:12) = eye(3)      
    M_Beam_hat(22:24,22:24) = eye(3)
    
! Calculation of the forces
    this%matrix12_object(i)%fqdyn_AB(:)  = 0.0d0
    this%matrix12_object(i)%fv_AB(:)     = 0.0d0
    !this%matrix12_object(i)%fqdyn_AB     = matmul(this%matrix12_object(i)%M_tilde,aAB)   
    !this%matrix12_object(i)%fv_AB        = matmul(this%matrix12_object(i)%M_tilde_hat,((qtn1AB-qtnAB)/deltat - 0.5d0*(vtn1AB+vtnAB)))- matmul(this%matrix12_object(i)%M_tilde, this%matrix12_object(i)%v_diss)
    this%matrix12_object(i)%fqdyn_AB      = matmul(this%matrix12_object(i)%M_Beam,aAB) 
    this%matrix12_object(i)%fv_AB         = matmul(M_Beam_hat,((qtn1AB-qtnAB)/deltat - 0.5d0*(vtn1AB+vtnAB)))
 
! Elements of the iteration matrix
    this%matrix12_object(i)%kqv(:,:)     = 0.0d0
    this%matrix12_object(i)%kvv(:,:)     = 0.0d0
    this%matrix12_object(i)%kvq(:,:)     = 0.0d0
    
    ! Rigid Body
    !this%matrix12_object(i)%kqv(:,:)     = this%matrix12_object(i)%M_tilde/deltat
    !this%matrix12_object(i)%kvv(:,:)     = -this%matrix12_object(i)%M_tilde_hat*0.5d0 - matmul(this%matrix12_object(i)%M_tilde,this%matrix12_object(i)%kvv_diss)
    !this%matrix12_object(i)%kvq(:,:)     = this%matrix12_object(i)%M_tilde_hat/deltat
    
    ! Beam 
    this%matrix12_object(i)%kqv(:,:)     =  this%matrix12_object(i)%M_Beam(:,:)/deltat
    this%matrix12_object(i)%kvv(:,:)     =  -0.5d0*M_Beam_hat(:,:)
    this%matrix12_object(i)%kvq(:,:)     =  M_Beam_hat(:,:)/deltat
           
            
end subroutine matrix12_mass
 
subroutine matrix12_damping(this, qtn1a, qtna, qtn1b, qtnb, vtn1a, vtn1b, vtna, vtnb, i)
   use my_math_structure, only: outer, eye, cross, skew
   use ieee_arithmetic

   implicit none
        
    class(matrix_12), intent(inout) :: this
    integer, intent(in)             :: i
    real(kind =8), dimension(12), intent(in)    :: qtn1a, qtna, qtn1b, qtnb, vtn1a, vtn1b, vtna, vtnb
    
    ! variables for damping functions from the beam element
    real(kind = 8), dimension(6)     :: strain_1, strain_2
    real(kind = 8), dimension(6,6)   :: C_in
    real(kind = 8)                   :: alpha_s, alpha_v, alpha_obj
    real(kind = 8), dimension(24)    :: v_1, v_2
    real(kind = 8), dimension(24,24) :: m
    real(kind = 8)                   :: normM_v_1, normM_v_2
    real(kind = 8)                   :: Dv, hv, Dv_v(24), hv_v(24) 
    
    ! variables for the modified damping function including the lie derivative
    real(kind = 8), dimension(6)     :: u1_dot, u2_dot, ua_dot, u_a 
    real(kind = 8), dimension(3)     :: av1a, av2a, av1b, av2b, ava          
    real(kind = 8), dimension(3)     :: d21a, d22a, d23a, phi2a          
    real(kind = 8), dimension(3)     :: d21b, d22b, d23b, phi2b         
    real(kind = 8), dimension(3)     :: d11a, d12a, d13a, phi1a         
    real(kind = 8), dimension(3)     :: d11b, d12b, d13b, phi1b   
    real(kind = 8), dimension(3)     :: d21a_dot, d22a_dot, d23a_dot, phi2a_dot          
    real(kind = 8), dimension(3)     :: d21b_dot, d22b_dot, d23b_dot, phi2b_dot         
    real(kind = 8), dimension(3)     :: d11a_dot, d12a_dot, d13a_dot, phi1a_dot         
    real(kind = 8), dimension(3)     :: d11b_dot, d12b_dot, d13b_dot, phi1b_dot   
    real(kind = 8), dimension(6,6)   :: Omega
    real(kind = 8)                   :: eta = 1.0d0
    
    ! Auxilliary variables to debug
    real(kind = 8), dimension(6)     :: Lie_u, diff, num, Omega_ua , res, f_diff   
    real(kind = 8)                   :: delta_t, D_analyt, Denum, Numerator

! ----------------------------------------------------------------------------------------
! 1. Dissipation from beam element
! ----------------------------------------------------------------------------------------
   
    strain_1 = this%matrix12_object(i)%strain1
    strain_2 = this%matrix12_object(i)%strain2
    u_a      = 0.5d0*(strain_2 + strain_1)
    
    alpha_s  = this%matrix12_property(i)%alpha_s
    alpha_v  = this%matrix12_property(i)%alpha_v
    
    C_in = this%matrix12_object(i)%c_in

!< Disspiation adopted from class_beam_element_dissipation -
    !< Stress dissipation: 
   if (alpha_s .ne. 0.0d0) then
      this%matrix12_object(i)%s_diss = 0.5d0*alpha_s*matmul(C_in,strain_2-strain_1)
      this%matrix12_object(i)%c_diss = 0.5d0*alpha_s*C_in
   end if
   
!< Dissipation adopted from class_beam_element_dissipation!
   !< Velocity dissipation: 
    if (alpha_v .ne. 0.0d0) then
        v_1(1:12)  = vtna
        v_1(13:24) = vtnb
        v_2(1:12)  = vtn1a
        v_2(13:24) = vtn1b
          
        m = this%matrix12_object(i)%M_tilde
        normM_v_1 = dsqrt(dot_product(v_1,matmul(m,v_1)))
        normM_v_2 = dsqrt(dot_product(v_2,matmul(m,v_2)))
        Dv = (normM_v_2 - normM_v_1);
        hv = (normM_v_2 + normM_v_1);
        if ((hv .eq. 0.0d0) .or. (normM_v_2 .eq. 0.0d0)) then
          this%matrix12_object(i)%v_diss(:) = 0.0d0
          this%matrix12_object(i)%kvv_diss  = 0.0d0
          return
        end if
        this%matrix12_object(i)%v_diss(:) = 0.5d0*alpha_v*(v_2+v_1)*Dv/hv
        Dv_v = matmul(m,v_2)/normM_v_2
        hv_v = Dv_v
        this%matrix12_object(i)%kvv_diss  = 0.5d0*alpha_v*Dv/hv*eye(24) + 0.5d0*alpha_v*outer(v_2+v_1,(Dv_v*hv - Dv*hv_v)/(hv*hv))  
    end if
    
    
!< ----------------------------------------------------------------------------------------
!< 2. Rayleigh damping from beam element built with lie derivatives 
!< ----------------------------------------------------------------------------------------
    
    if (this%matrix12_property(i)%alpha_obj .ne. 0.0d0) then 
       
        alpha_obj = this%matrix12_property(i)%alpha_obj 
        delta_t =  this%matrix12_object(i)%deltat
         
!< All positions and velocities to calculate the analyticaland numerical rates
        !< tn, node A
        phi1a  = qtna(1:3)
        d11a   = qtna(4:6)
        d12a   = qtna(7:9)
        d13a   = qtna(10:12)    
    
        !< tn+1, node A
        phi2a  = qtn1a(1:3)
        d21a   = qtn1a(4:6)
        d22a   = qtn1a(7:9)
        d23a   = qtn1a(10:12) 

        !< tn, node B
        phi1b  = qtnb(1:3)
        d11b   = qtnb(4:6)
        d12b   = qtnb(7:9)
        d13b   = qtnb(10:12)    
    
        !< tn+1, node B
        phi2b  = qtn1b(1:3)
        d21b   = qtn1b(4:6)
        d22b   = qtn1b(7:9)
        d23b   = qtn1b(10:12) 

        !< tn, node A
        phi1a_dot  = vtna(1:3)
        d11a_dot   = vtna(4:6)
        d12a_dot   = vtna(7:9)
        d13a_dot   = vtna(10:12)  
        
        !< tn+1, node A
        phi2a_dot  = vtn1a(1:3)
        d21a_dot   = vtn1a(4:6)
        d22a_dot   = vtn1a(7:9)
        d23a_dot   = vtn1a(10:12) 
 
        !< tn, node B
        phi1b_dot  = vtnb(1:3)
        d11b_dot   = vtnb(4:6)
        d12b_dot   = vtnb(7:9)
        d13b_dot   = vtnb(10:12)    
    
        !< tn+1, node B
        phi2b_dot  = vtn1b(1:3)
        d21b_dot   = vtn1b(4:6)
        d22b_dot   = vtn1b(7:9)
        d23b_dot   = vtn1b(10:12) 

        !< Analytical derivation of the deformation measure
        !< tn
        u1_dot(1) = dot_product(phi1b,(d11a_dot+d11b_dot)) + dot_product(phi1b_dot,(d11a+d11b)) - dot_product(phi1a,(d11a_dot+d11b_dot)) - dot_product(phi1a_dot,(d11a+d11b))
        u1_dot(2) = dot_product(phi1b,(d12a_dot+d12b_dot)) + dot_product(phi1b_dot,(d12a+d12b)) - dot_product(phi1a,(d12a_dot+d12b_dot)) - dot_product(phi1a_dot,(d12a+d12b))
        u1_dot(3) = dot_product(phi1b,(d13a_dot+d13b_dot)) + dot_product(phi1b_dot,(d13a+d13b)) - dot_product(phi1a,(d13a_dot+d13b_dot)) - dot_product(phi1a_dot,(d13a+d13b))
        u1_dot(4) = dot_product(d12a_dot,d13b) + dot_product(d13b_dot,d12a) - dot_product(d13a_dot,d12b) - dot_product(d12b_dot,d13a)
        u1_dot(5) = dot_product(d13a_dot,d11b) + dot_product(d11b_dot,d13a) - dot_product(d11a_dot,d13b) - dot_product(d13b_dot,d11a)
        u1_dot(6) = dot_product(d11a_dot,d12b) + dot_product(d12b_dot,d11a) - dot_product(d12a_dot,d11b) - dot_product(d11b_dot,d12a)
        u1_dot(1:6) = (1.0d0/(2.0d0*eta))*u1_dot(1:6)
        
        !< tn+1
        u2_dot(1) = dot_product(phi2b,(d21a_dot+d21b_dot)) + dot_product(phi2b_dot,(d21a+d21b)) - dot_product(phi2a,(d21a_dot+d21b_dot)) - dot_product(phi2a_dot,(d21a+d21b))
        u2_dot(2) = dot_product(phi2b,(d22a_dot+d22b_dot)) + dot_product(phi2b_dot,(d22a+d22b)) - dot_product(phi2a,(d22a_dot+d22b_dot)) - dot_product(phi2a_dot,(d22a+d22b))
        u2_dot(3) = dot_product(phi2b,(d23a_dot+d23b_dot)) + dot_product(phi2b_dot,(d23a+d23b)) - dot_product(phi2a,(d23a_dot+d23b_dot)) - dot_product(phi2a_dot,(d23a+d23b))
        u2_dot(4) = dot_product(d22a_dot,d23b) + dot_product(d23b_dot,d22a) - dot_product(d23a_dot,d22b) - dot_product(d22b_dot,d23a)
        u2_dot(5) = dot_product(d23a_dot,d21b) + dot_product(d21b_dot,d23a) - dot_product(d21a_dot,d23b) - dot_product(d23b_dot,d21a)
        u2_dot(6) = dot_product(d21a_dot,d22b) + dot_product(d22b_dot,d21a) - dot_product(d22a_dot,d21b) - dot_product(d21b_dot,d22a)
        u2_dot(1:6) = (1.0d0/(2.0d0*eta))*u2_dot(1:6)
                
        !< strain rates - averaged between tn+1 and tn
        ua_dot(1:6) = 0.5d0*(u1_dot(1:6)+u2_dot(1:6))
        
        !< calculation of the angular velocity
        !< tn
        !av1a(1) = dot_product(cross(d11a,d11a_dot) + cross(d12a,d12a_dot) + cross(d13a,d13a_dot), d11a) 
        !av1a(2) = dot_product(cross(d11a,d11a_dot) + cross(d12a,d12a_dot) + cross(d13a,d13a_dot), d12a) 
        !av1a(3) = dot_product(cross(d11a,d11a_dot) + cross(d12a,d12a_dot) + cross(d13a,d13a_dot), d13a) 
        !
        !av1b(1) = dot_product(cross(d11b,d11b_dot) + cross(d12b,d12b_dot) + cross(d13b,d13b_dot), d11b) 
        !av1b(2) = dot_product(cross(d11b,d11b_dot) + cross(d12b,d12b_dot) + cross(d13b,d13b_dot), d12b) 
        !av1b(3) = dot_product(cross(d11b,d11b_dot) + cross(d12b,d12b_dot) + cross(d13b,d13b_dot), d13b) 
        !
        !!< tn+1
        !av2a(1) = dot_product(cross(d21a,d21a_dot) + cross(d22a,d22a_dot) + cross(d23a,d23a_dot), d21a) 
        !av2a(2) = dot_product(cross(d21a,d21a_dot) + cross(d22a,d22a_dot) + cross(d23a,d23a_dot), d22a) 
        !av2a(3) = dot_product(cross(d21a,d21a_dot) + cross(d22a,d22a_dot) + cross(d23a,d23a_dot), d23a) 
        !
        !av2b(1) = dot_product(cross(d21b,d21b_dot) + cross(d22b,d22b_dot) + cross(d23b,d23b_dot), d21b) 
        !av2b(2) = dot_product(cross(d21b,d21b_dot) + cross(d22b,d22b_dot) + cross(d23b,d23b_dot), d22b) 
        !av2b(3) = dot_product(cross(d21b,d21b_dot) + cross(d22b,d22b_dot) + cross(d23b,d23b_dot), d23b) 
        
        !< angular velocities - averaged between tn+1 and tn 
        !ava(1:3) = 0.125d0*(av1a + av1b + av2a + av2b)
        !Omega(1:6,1:6) = 0.0d0
        !Omega(1:3,1:3) = skew(ava) 
        !Omega(4:6,4:6) = skew(ava) 
        !
        !!< Lie derivative
        !Lie_u = ua_dot - matmul(Omega,u_a)
        
        !<----------------------------------------------------------------
        !< Forces
            this%matrix12_object(i)%s_obj_diss = (0.5d0*alpha_obj) * matmul(C_in,delta_t*ua_dot)        
        !< Tangent
            this%matrix12_object(i)%c_obj_diss = (0.5d0*alpha_obj) * C_in
        !<----------------------------------------------------------------
        
        !< Forces neu
        !< Dissipationsfunktion einfach mit Ableitungen gebildet und in Zusammenhang aus Romero/Armero eingesetzt
        
            !D_analyt = sqrt(dot_product((u2_dot-u1_dot),matmul(C_in,(u2_dot-u1_dot)))) * sqrt(dot_product((u2_dot-u1_dot),matmul(C_in,(u2_dot-u1_dot))))
            !D_analyt = (0.5d0*alpha_obj*delta_t)*D_analyt            
            !Denum       = sqrt(dot_product((strain_2-strain_1),matmul(C_in,(strain_2-strain_1)))) * sqrt(dot_product((strain_2-strain_1),matmul(C_in,(strain_2-strain_1))))
            !
            !if (isnan(Denum) .or. Denum == 0.0d0) then
            !    this%matrix12_object(i)%s_obj_diss = 0.0d0
            !    this%matrix12_object(i)%c_obj_diss = 0.0d0
            !else
            !    this%matrix12_object(i)%s_obj_diss = (D_analyt/Denum)  * matmul(C_in,(strain_2-strain_1)) 
            !    this%matrix12_object(i)%c_obj_diss = (0.5d0*alpha_obj) * C_in
            !end if 
            
        
        !<----------------------------------------------------------------
        !< Comparing the lie derivative and the numerical derivative
            !this%matrix12_object(i)%num    = (strain_2 - strain_1)/delta_t
            !this%matrix12_object(i)%lie_u  = lie_u
            !this%matrix12_object(i)%ua_dot = ua_dot 
        !<----------------------------------------------------------------

    else
        this%matrix12_object(i)%s_obj_diss = 0.0d0
        this%matrix12_object(i)%c_obj_diss = 0.0d0
    end if
    



    return
end subroutine matrix12_damping




subroutine write_output(this, Lie, Num, diff, ua_dot,Omega_ua)
    implicit none
    class(matrix_12), intent(inout)            :: this
    real(kind = 8), dimension(6), intent(in)   :: Lie, Num, diff,ua_dot, Omega_ua
    character(len = 20)                        :: filename
    integer                     :: i 
    logical :: fileExist

    filename = "CE_output.txt"
    
    inquire(file=filename,exist=fileExist)
    if(fileExist) then
        open(unit=1234, file = filename, status = "replace")
    else
        open(unit=1234, file = filename, status = "new")
    end if
    do i = 1,6 
        write(1234,"(E15.8,2x,E15.8,2x,E15.8,2x,E15.8,2x,E15.8)") Lie(i), Num(i), diff(i), ua_dot(i), Omega_ua(i)
    end do
    close(1234)
    return 
end subroutine write_output

  

    end module class_matrix12    
