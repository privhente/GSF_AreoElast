module class_constraint_12
  
  use class_constraint
  use my_math_structure, only: eye, cross, skew, outer, kron, levi
  
  implicit none
    
  type, extends(constraint) :: constraint_12
    contains
      procedure :: gconstraint12  => gdummy12
      procedure :: dgqconstraint12 => dgqdummy12
      procedure :: dgvconstraint12 => dgvdummy12
      procedure :: kconstraint12  => kdummy12
  end type constraint_12
   
  type constraint_12_c   !! container type                                               
    class(constraint_12), pointer :: c
  end type
    
!> beam internal constraint class
  type, extends(constraint_12) :: constraint_12_internal
  contains
    procedure :: gconstraint12  => ginternal12
    procedure :: dgqconstraint12 => dgqinternal12
    procedure :: kconstraint12  => kinternal12
  end type constraint_12_internal
  
  interface constraint_12_internal
    module procedure create_constraint_12_internal
  end interface

!> beam reduced internal constraint class
  type, extends(constraint_12) :: constraint_12_internal_red
  contains
    procedure :: gconstraint12  => ginternal12_red
    procedure :: dgqconstraint12 => dgqinternal12_red
    procedure :: kconstraint12  => kinternal12_red
  end type constraint_12_internal_red
  
  interface constraint_12_internal_red
    module procedure create_constraint_12_internal_red
  end interface
  
!> beam spherical support constraint class
  type, extends(constraint_12) :: constraint_12_sphericalsupport
  contains
    procedure :: gconstraint12  => gsphericalsupport12
    procedure :: dgqconstraint12 => dgqsphericalsupport12
  end type constraint_12_sphericalsupport
  
  interface constraint_12_sphericalsupport
    module procedure create_constraint_12_sphericalsupport
  end interface
  
!> beam rigid support constraint class
  type, extends(constraint_12) :: constraint_12_rigidsupport
  contains
    procedure :: gconstraint12     => grigidsupport12
    procedure :: dgqconstraint12   => dgqrigidsupport12
  end type constraint_12_rigidsupport
  
  interface constraint_12_rigidsupport
    module procedure create_constraint_12_rigidsupport
  end interface
  
!> beam spherical joint constraint class
  type, extends(constraint_12) :: constraint_12_sphericaljoint
  contains
    procedure :: gconstraint12  => gsphericaljoint12
    procedure :: dgqconstraint12 => dgqsphericaljoint12
  end type constraint_12_sphericaljoint
  
  interface constraint_12_sphericaljoint
    module procedure create_constraint_12_sphericaljoint
  end interface
  
!> beam revolute joint constraint class
  type, extends(constraint_12) :: constraint_12_revolutejoint
  contains
    procedure :: gconstraint12  => grevolutejoint12
    procedure :: dgqconstraint12 => dgqrevolutejoint12
    procedure :: kconstraint12  => krevolutejoint12
  end type constraint_12_revolutejoint
  
  interface constraint_12_revolutejoint
    module procedure create_constraint_12_revolutejoint
  end interface
  
!> beam rotor revolute joint constraint class
  type, extends(constraint_12) :: constraint_12_rotorrevolutejoint
  contains
    procedure :: gconstraint12  => grotorrevolutejoint12
    procedure :: dgqconstraint12 => dgqrotorrevolutejoint12
    procedure :: kconstraint12  => krotorrevolutejoint12
  end type constraint_12_rotorrevolutejoint
  
  interface constraint_12_rotorrevolutejoint
    module procedure create_constraint_12_rotorrevolutejoint
  end interface
  
!> beam new revolute joint constraint class
  type, extends(constraint_12) :: constraint_12_revolutejoint_2
  contains
    procedure :: gconstraint12  => grevolutejoint12_2
    procedure :: dgqconstraint12 => dgqrevolutejoint12_2
    procedure :: kconstraint12  => krevolutejoint12_2
  end type constraint_12_revolutejoint_2
  
  interface constraint_12_revolutejoint_2
    module procedure create_constraint_12_revolutejoint_2
  end interface
  
!> beam new revolute joint constraint class
  type, extends(constraint_12) :: constraint_12_revolutejoint_3
  contains
    procedure :: gconstraint12  => grevolutejoint12_3
    procedure :: dgqconstraint12 => dgqrevolutejoint12_3
    procedure :: kconstraint12  => krevolutejoint12_3
  end type constraint_12_revolutejoint_3
  
  interface constraint_12_revolutejoint_3
    module procedure create_constraint_12_revolutejoint_3
  end interface


!> beam rigid connection constraint class
  type, extends(constraint_12) :: constraint_12_rigidconnection
  contains
    procedure :: gconstraint12   => grigidconnection12
    procedure :: dgqconstraint12 => dgqrigidconnection12
    procedure :: kconstraint12   => krigidconnection12
  end type constraint_12_rigidconnection
  
  interface constraint_12_rigidconnection
    module procedure create_constraint_12_rigidconnection
  end interface
  
!> beam masspoint connection constraint class
  type, extends(constraint_12) :: constraint_12_masspointconnection
  contains
    procedure :: gconstraint12  => gmasspointconnection12
    procedure :: dgqconstraint12 => dgqmasspointconnection12
  end type constraint_12_masspointconnection
  
  interface constraint_12_masspointconnection
    module procedure create_constraint_12_masspointconnection
  end interface
  
!> beam simple support constraint class
  type, extends(constraint_12) :: constraint_12_simplesupport_global
  contains
    procedure :: gconstraint12  => gsimplesupport_global
    procedure :: dgqconstraint12 => dgqsimplesupport_global
  end type constraint_12_simplesupport_global
  
  interface constraint_12_simplesupport_global
    module procedure create_constraint_12_simplesupport_global
  end interface  

!> beam simple support constraint class
  type, extends(constraint_12) :: constraint_12_simplesupport_local
  contains
    procedure :: gconstraint12     => gsimplesupport_local
    procedure :: dgqconstraint12   => dgqsimplesupport_local
    procedure :: kconstraint12     => ksimplesupport_local
  end type constraint_12_simplesupport_local
  
  interface constraint_12_simplesupport_local
    module procedure create_constraint_12_simplesupport_local
  end interface  

!> Node 12 angular velocity around a defined axis constraint class
  type, extends(constraint_12) :: constraint_12_angularvelocity_global
  contains
    procedure :: gconstraint12  => gangularvelocity_global
    procedure :: dgqconstraint12 => dgqangularvelocity_global
    procedure :: dgvconstraint12 => dgvangularvelocity_global
    procedure :: kconstraint12   => kangularvelocity_global
  end type constraint_12_angularvelocity_global
  
  interface constraint_12_angularvelocity_global
    module procedure create_constraint_12_angularvelocity_global
  end interface
  
!> Node 12 angular velocity around a defined axis constraint class
  type, extends(constraint_12) :: constraint_12_angularvelocity_local
  contains
    procedure :: gconstraint12  => gangularvelocity_local
    procedure :: dgqconstraint12 => dgqangularvelocity_local
    procedure :: dgvconstraint12 => dgvangularvelocity_local
    procedure :: kconstraint12   => kangularvelocity_local
  end type constraint_12_angularvelocity_local
  
  interface constraint_12_angularvelocity_local
    module procedure create_constraint_12_angularvelocity_local
  end interface

!> Node 12 rotation around a defined axis constraint class
  type, extends(constraint_12) :: constraint_12_rotation_local
  contains
    procedure :: gconstraint12  => grotation_local
    procedure :: dgqconstraint12 => dgqrotation_local
    procedure :: kconstraint12   => krotation_local
  end type constraint_12_rotation_local
  
  interface constraint_12_rotation_local
    module procedure create_constraint_12_rotation_local
  end interface

  !> Node 12 rotation around a defined axis constraint class relative to another node
  type, extends(constraint_12) :: constraint_12_relative_rotation_local
  contains
    procedure :: gconstraint12  => grelative_rotation_local
    procedure :: dgqconstraint12 => dgqrelative_rotation_local
    procedure :: kconstraint12 => krelative_rotation_local
  end type constraint_12_relative_rotation_local
  
  interface constraint_12_relative_rotation_local
    module procedure create_constraint_12_relative_rotation_local
  end interface
  
  
!> Node 12 rotation around a defined global axis 
  type, extends(constraint_12) :: constraint_12_rotation_global
  contains
    procedure :: gconstraint12  => grotation_global
    procedure :: dgqconstraint12 => dgqrotation_global
  end type constraint_12_rotation_global
  
  interface constraint_12_rotation_global
    module procedure create_constraint_12_rotation_global
  end interface

! SUBROUTINES  
  contains

! ====================================================================================================================
! DUMMY CONSTRAINT  
! ====================================================================================================================
  subroutine gdummy12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    g(:) = 0.0d0
  end subroutine gdummy12

  subroutine dgqdummy12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    dgq(:,:) = 0.0d0
  end subroutine dgqdummy12

  subroutine dgvdummy12(this, q1_t, q2_t, q1_0, q2_0, dgv, qopt1, qopt2)
    implicit none
    class(constraint_12), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgv(:,:)
    dgv(:,:) = 0.0d0
  end subroutine dgvdummy12  

  subroutine kdummy12(this, lambda, k, qopt)
    implicit none
    class(constraint_12), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    k(:,:) = 0.0d0
  end subroutine kdummy12

! ====================================================================================================================
!> Defining constraint variables
! ====================================================================================================================
  type(constraint_12_internal) function create_constraint_12_internal()
    create_constraint_12_internal%rankcount = 6
    create_constraint_12_internal%coordinates = 1
    create_constraint_12_internal%stiffness = .TRUE.
    create_constraint_12_internal%boolean_q  = .TRUE.
    create_constraint_12_internal%boolean_v  = .FALSE.
    return
  end function
  
  type(constraint_12_internal_red) function create_constraint_12_internal_red()
    create_constraint_12_internal_red%rankcount = 6
    create_constraint_12_internal_red%coordinates = 1
    create_constraint_12_internal_red%stiffness = .TRUE.
    create_constraint_12_internal_red%boolean_q  = .TRUE.
    create_constraint_12_internal_red%boolean_v  = .FALSE.
    return
  end function

  type(constraint_12_sphericalsupport) function create_constraint_12_sphericalsupport()
    create_constraint_12_sphericalsupport%rankcount = 3
    create_constraint_12_sphericalsupport%coordinates = 1
    create_constraint_12_sphericalsupport%stiffness = .FALSE.
    create_constraint_12_sphericalsupport%boolean_q  = .TRUE.
    create_constraint_12_sphericalsupport%boolean_v  = .FALSE.
    return
  end function
  
  type(constraint_12_rigidsupport) function create_constraint_12_rigidsupport()
    create_constraint_12_rigidsupport%rankcount = 12
    create_constraint_12_rigidsupport%coordinates = 1
    create_constraint_12_rigidsupport%stiffness = .FALSE.
    create_constraint_12_rigidsupport%boolean_q  = .TRUE.
    create_constraint_12_rigidsupport%boolean_v  = .FALSE.
    return
  end function
  
  type(constraint_12_simplesupport_global) function create_constraint_12_simplesupport_global()
    create_constraint_12_simplesupport_global%rankcount = 1
    create_constraint_12_simplesupport_global%coordinates = 1
    create_constraint_12_simplesupport_global%stiffness = .FALSE.
    create_constraint_12_simplesupport_global%boolean_q  = .TRUE.
    create_constraint_12_simplesupport_global%boolean_v  = .FALSE.
    return
  end function
  
  type(constraint_12_simplesupport_local) function create_constraint_12_simplesupport_local()
    create_constraint_12_simplesupport_local%rankcount = 1
    create_constraint_12_simplesupport_local%coordinates = 1
    create_constraint_12_simplesupport_local%stiffness = .TRUE.
    create_constraint_12_simplesupport_local%boolean_q  = .TRUE.
    create_constraint_12_simplesupport_local%boolean_v  = .FALSE.
    return
  end function
  
  type(constraint_12_sphericaljoint) function create_constraint_12_sphericaljoint()
    create_constraint_12_sphericaljoint%rankcount = 3
    create_constraint_12_sphericaljoint%coordinates = 2
    create_constraint_12_sphericaljoint%stiffness = .FALSE.
    create_constraint_12_sphericaljoint%boolean_q  = .TRUE.
    create_constraint_12_sphericaljoint%boolean_v  = .FALSE.
    return
  end function
  
  type(constraint_12_revolutejoint) function create_constraint_12_revolutejoint()
    create_constraint_12_revolutejoint%rankcount = 5
    create_constraint_12_revolutejoint%coordinates = 2
    create_constraint_12_revolutejoint%stiffness = .TRUE.
    create_constraint_12_revolutejoint%boolean_q  = .TRUE.
    create_constraint_12_revolutejoint%boolean_v  = .FALSE.
    return
  end function

  type(constraint_12_rotorrevolutejoint) function create_constraint_12_rotorrevolutejoint()
    create_constraint_12_rotorrevolutejoint%rankcount = 6
    create_constraint_12_rotorrevolutejoint%coordinates = 2
    create_constraint_12_rotorrevolutejoint%stiffness = .TRUE.
    create_constraint_12_rotorrevolutejoint%boolean_q  = .TRUE.
    create_constraint_12_rotorrevolutejoint%boolean_v  = .FALSE.
    return
  end function

  type(constraint_12_revolutejoint_2) function create_constraint_12_revolutejoint_2()
    create_constraint_12_revolutejoint_2%rankcount = 6
    create_constraint_12_revolutejoint_2%coordinates = 2
    create_constraint_12_revolutejoint_2%stiffness = .TRUE.
    create_constraint_12_revolutejoint_2%boolean_q  = .TRUE.
    create_constraint_12_revolutejoint_2%boolean_v  = .FALSE.
    return
    end function
    
  type(constraint_12_revolutejoint_3) function create_constraint_12_revolutejoint_3()
    create_constraint_12_revolutejoint_3%rankcount = 5
    create_constraint_12_revolutejoint_3%coordinates = 2
    create_constraint_12_revolutejoint_3%stiffness = .TRUE.
    create_constraint_12_revolutejoint_3%boolean_q  = .TRUE.
    create_constraint_12_revolutejoint_3%boolean_v  = .FALSE.
    return
  end function

  type(constraint_12_rigidconnection) function create_constraint_12_rigidconnection()
    create_constraint_12_rigidconnection%rankcount = 12
    create_constraint_12_rigidconnection%coordinates = 2
    create_constraint_12_rigidconnection%stiffness = .TRUE.
    create_constraint_12_rigidconnection%boolean_q  = .TRUE.
    create_constraint_12_rigidconnection%boolean_v  = .FALSE.
    return
  end function
  
  type(constraint_12_masspointconnection) function create_constraint_12_masspointconnection()
    create_constraint_12_masspointconnection%rankcount = 12
    create_constraint_12_masspointconnection%coordinates = 2
    create_constraint_12_masspointconnection%stiffness = .FALSE.
    create_constraint_12_masspointconnection%boolean_q  = .TRUE.
    create_constraint_12_masspointconnection%boolean_v  = .FALSE.
    return
  end function

  type(constraint_12_angularvelocity_global) function create_constraint_12_angularvelocity_global()
    create_constraint_12_angularvelocity_global%rankcount = 1
    create_constraint_12_angularvelocity_global%coordinates = 1
    create_constraint_12_angularvelocity_global%stiffness = .TRUE.
    create_constraint_12_angularvelocity_global%boolean_q  = .TRUE.
    create_constraint_12_angularvelocity_global%boolean_v  = .TRUE.
    return
  end function
  
  type(constraint_12_angularvelocity_local) function create_constraint_12_angularvelocity_local()
    create_constraint_12_angularvelocity_local%rankcount = 1
    create_constraint_12_angularvelocity_local%coordinates = 1
    create_constraint_12_angularvelocity_local%stiffness = .TRUE.
    create_constraint_12_angularvelocity_local%boolean_q  = .TRUE.
    create_constraint_12_angularvelocity_local%boolean_v  = .TRUE.
    return
  end function

  type(constraint_12_rotation_local) function create_constraint_12_rotation_local()
    create_constraint_12_rotation_local%rankcount = 1
    create_constraint_12_rotation_local%coordinates = 1
    create_constraint_12_rotation_local%stiffness = .TRUE.
    create_constraint_12_rotation_local%boolean_q  = .TRUE.
    create_constraint_12_rotation_local%boolean_v  = .FALSE.
    return
    end function
    
  type(constraint_12_relative_rotation_local) function create_constraint_12_relative_rotation_local()
    create_constraint_12_relative_rotation_local%rankcount = 12
    create_constraint_12_relative_rotation_local%coordinates = 2
    create_constraint_12_relative_rotation_local%stiffness = .TRUE.
    create_constraint_12_relative_rotation_local%boolean_q  = .TRUE.
    create_constraint_12_relative_rotation_local%boolean_v  = .FALSE.
    return
  end function

  type(constraint_12_rotation_global) function create_constraint_12_rotation_global()
    create_constraint_12_rotation_global%rankcount = 1
    create_constraint_12_rotation_global%coordinates = 1
    create_constraint_12_rotation_global%stiffness = .FALSE.
    create_constraint_12_rotation_global%boolean_q  = .TRUE.
    create_constraint_12_rotation_global%boolean_v  = .FALSE.
    return
  end function

! ====================================================================================================================
! INTERNAL CONSTRAINT - LENGTH AND ORTHOGONALITY  
! ====================================================================================================================
  subroutine ginternal12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none

    class(constraint_12_internal), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(3) :: d1_t, d2_t, d3_t
    real(kind = 8), dimension(3) :: d1_0, d2_0, d3_0
    
    d1_t = q1_t( 4: 6)
    d2_t = q1_t( 7: 9)
    d3_t = q1_t(10:12)
  
    d1_0 = q1_0( 4: 6)
    d2_0 = q1_0( 7: 9)
    d3_0 = q1_0(10:12)
  
    g(1) = 0.5d0*(dot_product(d1_t, d1_t)-dot_product(d1_0, d1_0))
    g(2) = 0.5d0*(dot_product(d2_t, d2_t)-dot_product(d2_0, d2_0))
    g(3) = 0.5d0*(dot_product(d3_t, d3_t)-dot_product(d3_0, d3_0))
    g(4) = dot_product(d1_t, d2_t)-dot_product(d1_0, d2_0)
    g(5) = dot_product(d1_t, d3_t)-dot_product(d1_0, d3_0)
    g(6) = dot_product(d2_t, d3_t)-dot_product(d2_0, d3_0)
    
  end subroutine ginternal12

  subroutine dgqinternal12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    
    implicit none

    class(constraint_12_internal), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(3) :: d1_t, d2_t, d3_t

    d1_t = q1_t( 4: 6)
    d2_t = q1_t( 7: 9)
    d3_t = q1_t(10:12)
   
    dgq( :, :) = 0.0d0

    dgq( 1, 4: 6) = d1_t
    dgq( 2, 7: 9) = d2_t
    dgq( 3,10:12) = d3_t

    dgq( 4, 4: 6) = d2_t
    dgq( 4, 7: 9) = d1_t

    dgq( 5, 4: 6) = d3_t
    dgq( 5,10:12) = d1_t

    dgq( 6, 7: 9) = d3_t
    dgq( 6,10:12) = d2_t
    
  end subroutine dgqinternal12

  subroutine kinternal12(this, lambda, k, qopt)

    implicit none

    class(constraint_12_internal), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    k(:, :) = 0.0d0

    k( 4: 6, 4: 6) = lambda(1)*i_33
    k( 7: 9, 7: 9) = lambda(2)*i_33
    k(10:12,10:12) = lambda(3)*i_33

    k( 4: 6, 7: 9) = lambda(4)*i_33
    k( 7: 9, 4: 6) = lambda(4)*i_33

    k( 4: 6,10:12) = lambda(5)*i_33
    k(10:12, 4: 6) = lambda(5)*i_33

    k( 7: 9,10:12) = lambda(6)*i_33
    k(10:12, 7: 9) = lambda(6)*i_33
    
  end subroutine kinternal12

! ====================================================================================================================
! REDUCED INTERNAL CONSTRAINT
! ====================================================================================================================
  subroutine ginternal12_red(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none

    class(constraint_12_internal_red), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(3) :: d1_t, d2_t, d3_t
    real(kind = 8), dimension(3) :: d1_0, d2_0, d3_0
    
    d1_t = q1_t( 4: 6)
    d2_t = q1_t( 7: 9)
    d3_t = q1_t(10:12)
  
    d1_0 = q1_0( 4: 6)
    d2_0 = q1_0( 7: 9)
    d3_0 = q1_0(10:12)
  
    g(1) = 0.5d0*(dot_product(d1_t, d1_t)-dot_product(d1_0, d1_0))
    g(2) = 0.5d0*(dot_product(d2_t, d2_t)-dot_product(d2_0, d2_0))
    g(3) = 0.5d0*(dot_product(d3_t, d3_t)-dot_product(d3_0, d3_0))
    g(4) = dot_product(d1_t, d2_t)-dot_product(d1_0, d2_0)
    g(5) = dot_product(d1_t, d3_t)-dot_product(d1_0, d3_0)
    g(6) = dot_product(d2_t, d3_t)-dot_product(d2_0, d3_0)
    
  end subroutine ginternal12_red

  subroutine dgqinternal12_red(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    
    implicit none

    class(constraint_12_internal_red), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(3) :: d1_t, d2_t, d3_t

    d1_t = q1_t( 4: 6)
    d2_t = q1_t( 7: 9)
    d3_t = q1_t(10:12)
   
    dgq( :, :) = 0.0d0

    dgq( 1, 4: 6) = d1_t
    dgq( 2, 7: 9) = d2_t
    dgq( 3,10:12) = d3_t

    dgq( 4, 4: 6) = d2_t
    dgq( 4, 7: 9) = d1_t

    dgq( 5, 4: 6) = d3_t
    dgq( 5,10:12) = d1_t

    dgq( 6, 7: 9) = d3_t
    dgq( 6,10:12) = d2_t
    
  end subroutine dgqinternal12_red

  subroutine kinternal12_red(this, lambda, k, qopt)

    implicit none

    class(constraint_12_internal_red), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    k(:, :) = 0.0d0

    !k( 4: 6, 4: 6) = lambda(1)*i_33
    !k( 7: 9, 7: 9) = lambda(2)*i_33
    !k(10:12,10:12) = lambda(3)*i_33

    !k( 4: 6, 7: 9) = lambda(4)*i_33
    !k( 7: 9, 4: 6) = lambda(4)*i_33
    !
    !k( 4: 6,10:12) = lambda(5)*i_33
    !k(10:12, 4: 6) = lambda(5)*i_33
    !
    !k( 7: 9,10:12) = lambda(6)*i_33
    !k(10:12, 7: 9) = lambda(6)*i_33
    !
  end subroutine kinternal12_red
  
! ====================================================================================================================
! SPHERICAL SUPPORT
! ====================================================================================================================
  subroutine gsphericalsupport12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)

    implicit none

    class(constraint_12_sphericalsupport), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    
    g(1:3) = (q1_t(1:3)+this%phi(1, 1)*q1_t(4:6)+this%phi(2, 1)*q1_t(7:9)+this%phi(3, 1)*q1_t(10:12))&
            -(q1_0(1:3)+this%phi(1, 1)*q1_0(4:6)+this%phi(2, 1)*q1_0(7:9)+this%phi(3, 1)*q1_0(10:12))
    
  end subroutine gsphericalsupport12

  subroutine dgqsphericalsupport12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)

    implicit none

    class(constraint_12_sphericalsupport), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    
    dgq(:, :) = 0.d0

    dgq(1:3,  1: 3) = i_33
    dgq(1:3,  4: 6) = i_33*this%phi(1, 1)
    dgq(1:3,  7: 9) = i_33*this%phi(2, 1)
    dgq(1:3, 10:12) = i_33*this%phi(3, 1)
    
  end subroutine dgqsphericalsupport12

! ====================================================================================================================
! RIGID SUPPORT
! ====================================================================================================================
  subroutine grigidsupport12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)

  implicit none

    class(constraint_12_rigidsupport), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)

    g(1:3) = (q1_t(1:3)+this%phi(1, 1)*q1_t(4:6)+this%phi(2, 1)*q1_t(7:9)+this%phi(3, 1)*q1_t(10:12))&
            -(q1_0(1:3)+this%phi(1, 1)*q1_0(4:6)+this%phi(2, 1)*q1_0(7:9)+this%phi(3, 1)*q1_0(10:12))
    
    g( 4:12) = q1_t(4:12)-q1_0(4:12)

  end subroutine grigidsupport12

    subroutine dgqrigidsupport12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)

    implicit none

    class(constraint_12_rigidsupport), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    
    dgq(:,:) = 0.0d0
    
    dgq(1:3,  1: 3) = i_33
    dgq(1:3,  4: 6) = i_33*this%phi(1, 1)
    dgq(1:3,  7: 9) = i_33*this%phi(2, 1)
    dgq(1:3, 10:12) = i_33*this%phi(3, 1)

    dgq( 4: 6,  4: 6) = i_33
    dgq( 7: 9,  7: 9) = i_33
    dgq(10:12, 10:12) = i_33

  end subroutine dgqrigidsupport12

! ====================================================================================================================
! SIMPLESUPPORT global
! ====================================================================================================================
    subroutine gsimplesupport_global(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)

    implicit none

    class(constraint_12_simplesupport_global), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8) :: u_n(3), u_tilde(3), x_n(3), x_0(3)
    
    x_n = q1_t(1:3) + this%phi(1,1)*q1_t(4:6) + this%phi(2,1)*q1_t(7:9) + this%phi(3,1)*q1_t(10:12)
    x_0 = q1_0(1:3) + this%phi(1,1)*q1_0(4:6) + this%phi(2,1)*q1_0(7:9) + this%phi(3,1)*q1_0(10:12)
    u_n = x_n - x_0
    
    u_tilde = this%amplitude_translation*this%dir
    
    g = dot_product(u_n - u_tilde, this%dir)
    
  end subroutine gsimplesupport_global

  subroutine dgqsimplesupport_global(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)

    implicit none

    class(constraint_12_simplesupport_global), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8) :: A_p(3,12)
    
    A_p(:,:)       = 0.0d0
    A_p(1:3,1:3)   = i_33
    A_p(1:3,4:6)   = i_33*this%phi(1,1)
    A_p(1:3,7:9)   = i_33*this%phi(2,1)
    A_p(1:3,10:12) = i_33*this%phi(3,1)
    
    dgq(1,:) = 0.0d0
    dgq(1,1:12) = matmul(this%dir,A_p)
    
  end subroutine dgqsimplesupport_global

! ====================================================================================================================
! SIMPLESUPPORT local
! ====================================================================================================================
    subroutine gsimplesupport_local(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)

    implicit none

    class(constraint_12_simplesupport_local), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8) :: u_n(3), u_tilde(3), x_n(3), x_0(3), n(3)
    
    x_n = q1_t(1:3) + this%phi(1,1)*q1_t(4:6) + this%phi(2,1)*q1_t(7:9) + this%phi(3,1)*q1_t(10:12)
    x_0 = q1_0(1:3) + this%phi(1,1)*q1_0(4:6) + this%phi(2,1)*q1_0(7:9) + this%phi(3,1)*q1_0(10:12)
    u_n = x_n - x_0
    n   = this%dir(1)*q1_t(4:6) + this%dir(2)*q1_t(7:9) + this%dir(3)*q1_t(10:12)
    
    u_tilde = this%amplitude_translation*n
    
    g = dot_product(u_n - u_tilde, n)
    
  end subroutine gsimplesupport_local

  subroutine dgqsimplesupport_local(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)

    implicit none

    class(constraint_12_simplesupport_local), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8) :: A_p(3,12), A_n(3,12), u_n(3), u_tilde(3), x_n(3), x_0(3), n(3)
    
    x_n = q1_t(1:3) + this%phi(1,1)*q1_t(4:6) + this%phi(2,1)*q1_t(7:9) + this%phi(3,1)*q1_t(10:12)
    x_0 = q1_0(1:3) + this%phi(1,1)*q1_0(4:6) + this%phi(2,1)*q1_0(7:9) + this%phi(3,1)*q1_0(10:12)
    
    n       = this%dir(1)*q1_t(4:6) + this%dir(2)*q1_t(7:9) + this%dir(3)*q1_t(10:12)
    u_n     = x_n - x_0
    
    A_p(:,:)       = 0.0d0
    A_p(1:3,1:3)   = i_33
    A_p(1:3,4:6)   = i_33*this%phi(1,1)
    A_p(1:3,7:9)   = i_33*this%phi(2,1)
    A_p(1:3,10:12) = i_33*this%phi(3,1)
    
    A_n(:,:)       = 0.0d0
    A_n(1:3,1:3)   = i_33*0.0d0
    A_n(1:3,4:6)   = i_33*this%dir(1)
    A_n(1:3,7:9)   = i_33*this%dir(2)
    A_n(1:3,10:12) = i_33*this%dir(3)
    
    dgq(1,:) = 0.0d0
    dgq(1,1:12) = matmul(n,A_p) + matmul(u_n,A_n)
    
  end subroutine dgqsimplesupport_local

  subroutine ksimplesupport_local(this, lambda, k, qopt)

    implicit none

    class(constraint_12_simplesupport_local), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    real(kind = 8) :: A_p(3,12), A_n(3,12), n(3)
    
    A_p(:,:)       = 0.0d0
    A_p(1:3,1:3)   = i_33
    A_p(1:3,4:6)   = i_33*this%phi(1,1)
    A_p(1:3,7:9)   = i_33*this%phi(2,1)
    A_p(1:3,10:12) = i_33*this%phi(3,1)
    
    A_n(:,:)       = 0.0d0
    A_n(1:3,1:3)   = i_33*0.0d0
    A_n(1:3,4:6)   = i_33*this%dir(1)
    A_n(1:3,7:9)   = i_33*this%dir(2)
    A_n(1:3,10:12) = i_33*this%dir(3)   
    
    k(:, :)      = 0.0d0
    k(1:12,1:12) = lambda(1)*( matmul(transpose(A_p),A_n) + matmul(transpose(A_n),A_p) )
    
  end subroutine ksimplesupport_local
  
! ====================================================================================================================
! SPHERICAL JOINT
! ====================================================================================================================  
  subroutine gsphericaljoint12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)

    implicit none

    class(constraint_12_sphericaljoint), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(3) :: phi10, phi20, phi1t, phi2t
    
    phi10  = q1_0(1:3) + this%phi(1,1)*q1_0(4:6) + this%phi(2,1)*q1_0(7:9) + this%phi(3,1)*q1_0(10:12)
    phi20  = q2_0(1:3) + this%phi(1,2)*q2_0(4:6) + this%phi(2,2)*q2_0(7:9) + this%phi(3,2)*q2_0(10:12)
    phi1t  = q1_t(1:3) + this%phi(1,1)*q1_t(4:6) + this%phi(2,1)*q1_t(7:9) + this%phi(3,1)*q1_t(10:12)
    phi2t  = q2_t(1:3) + this%phi(1,2)*q2_t(4:6) + this%phi(2,2)*q2_t(7:9) + this%phi(3,2)*q2_t(10:12)

    g(:)   = 0.0d0
    g(1:3) = (phi2t - phi1t) - (phi20 - phi10)
    
  end subroutine gsphericaljoint12

  subroutine dgqsphericaljoint12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)

    implicit none

    class(constraint_12_sphericaljoint), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    
    dgq(:,:)        = 0.0d0
    dgq(1:3,  1: 3) = -i_33
    dgq(1:3,  4: 6) = -i_33*this%phi(1, 1)
    dgq(1:3,  7: 9) = -i_33*this%phi(2, 1)
    dgq(1:3, 10:12) = -i_33*this%phi(3, 1)

    dgq(1:3, 13:15) = i_33
    dgq(1:3, 16:18) = i_33*this%phi(1, 2)
    dgq(1:3, 19:21) = i_33*this%phi(2, 2)
    dgq(1:3, 22:24) = i_33*this%phi(3, 2)
    
  end subroutine dgqsphericaljoint12

! ====================================================================================================================
! REVOLUTE JOINT
! ====================================================================================================================
  subroutine grevolutejoint12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_revolutejoint), intent(inout) :: this
    !in allen subroutinen g... auf inout \E4ndern, constraint 12
    ! weil \FCber dummy routine muss in allen g.. gleich sein
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) :: ns, ns0
    
    !< roation axis defined related to node 1 (dir in director 1 cos)
    ns  =   (q2_t(1:3)+this%phi(1, 2)*q2_t(4:6)+this%phi(2, 2)*q2_t(7:9)+this%phi(3, 2)*q2_t(10:12))& 
          - (q1_t(1:3)+this%phi(1, 1)*q1_t(4:6)+this%phi(2, 1)*q1_t(7:9)+this%phi(3, 1)*q1_t(10:12))
    
    ns0 =   (q2_0(1:3)+this%phi(1, 2)*q2_0(4:6)+this%phi(2, 2)*q2_0(7:9)+this%phi(3, 2)*q2_0(10:12))& 
          - (q1_0(1:3)+this%phi(1, 1)*q1_0(4:6)+this%phi(2, 1)*q1_0(7:9)+this%phi(3, 1)*q1_0(10:12))

    this%dir_gl = ns/norm2(ns) ! in global KOS to store for later use, if necessary
    
    ! conservation of angle and distance between distance vector at t and t0
    g(1) =  dot_product(ns,q1_t(4:6))   - dot_product(ns0,q1_0(4:6))
    g(2) =  dot_product(ns,q1_t(7:9))   - dot_product(ns0,q1_0(7:9))
    g(3) =  dot_product(ns,q1_t(10:12)) - dot_product(ns0,q1_0(10:12))
    
    ! conservation of angle between axis and directors of node 2
    g(4) = dot_product(ns, q2_t( 4: 6))-dot_product(ns0, q2_0( 4: 6))
    g(5) = dot_product(ns, q2_t(10:12))-dot_product(ns0, q2_0(10:12))

  end subroutine grevolutejoint12

  subroutine dgqrevolutejoint12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_revolutejoint), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(3) :: ns
    
    ns  =   (q2_t(1:3)+this%phi(1, 2)*q2_t(4:6)+this%phi(2, 2)*q2_t(7:9)+this%phi(3, 2)*q2_t(10:12)) & 
          - (q1_t(1:3)+this%phi(1, 1)*q1_t(4:6)+this%phi(2, 1)*q1_t(7:9)+this%phi(3, 1)*q1_t(10:12))
            
    dgq(:, :)     = 0.0d0
    
    dgq(1,  1: 3) =     - q1_t(4:6)
    dgq(1,  4: 6) =  ns - q1_t(4:6)*this%phi(1, 1)
    dgq(1,  7: 9) =     - q1_t(4:6)*this%phi(2, 1)
    dgq(1, 10:12) =     - q1_t(4:6)*this%phi(3, 1)
    dgq(1, 13:15) =       q1_t(4:6)
    dgq(1, 16:18) =       q1_t(4:6)*this%phi(1, 2)
    dgq(1, 19:21) =       q1_t(4:6)*this%phi(2, 2)
    dgq(1, 22:24) =       q1_t(4:6)*this%phi(3, 2)
    
    dgq(2,  1: 3) =     - q1_t(7:9)
    dgq(2,  4: 6) =     - q1_t(7:9)*this%phi(1, 1)
    dgq(2,  7: 9) =  ns - q1_t(7:9)*this%phi(2, 1)
    dgq(2, 10:12) =     - q1_t(7:9)*this%phi(3, 1)
    dgq(2, 13:15) =       q1_t(7:9)
    dgq(2, 16:18) =       q1_t(7:9)*this%phi(1, 2)
    dgq(2, 19:21) =       q1_t(7:9)*this%phi(2, 2) 
    dgq(2, 22:24) =       q1_t(7:9)*this%phi(3, 2)
    
    dgq(3,  1: 3) =     - q1_t(10:12)
    dgq(3,  4: 6) =     - q1_t(10:12)*this%phi(1, 1)
    dgq(3,  7: 9) =     - q1_t(10:12)*this%phi(2, 1)
    dgq(3, 10:12) =  ns - q1_t(10:12)*this%phi(3, 1) 
    dgq(3, 13:15) =       q1_t(10:12)
    dgq(3, 16:18) =       q1_t(10:12)*this%phi(1, 2)
    dgq(3, 19:21) =       q1_t(10:12)*this%phi(2, 2) 
    dgq(3, 22:24) =       q1_t(10:12)*this%phi(3, 2)
!    
    dgq(4,  1: 3) =     - q2_t(4:6)
    dgq(4,  4: 6) =     - q2_t(4:6)*this%phi(1, 1)
    dgq(4,  7: 9) =     - q2_t(4:6)*this%phi(2, 1)
    dgq(4, 10:12) =     - q2_t(4:6)*this%phi(3, 1)
    dgq(4, 13:15) =       q2_t(4:6)
    dgq(4, 16:18) =  ns + q2_t(4:6)*this%phi(1, 2)
    dgq(4, 19:21) =       q2_t(4:6)*this%phi(2, 2)
    dgq(4, 22:24) =       q2_t(4:6)*this%phi(3, 2)
    
    dgq(5,  1: 3) =     - q2_t(10:12)
    dgq(5,  4: 6) =     - q2_t(10:12)*this%phi(1, 1)
    dgq(5,  7: 9) =     - q2_t(10:12)*this%phi(2, 1)
    dgq(5, 10:12) =     - q2_t(10:12)*this%phi(3, 1)
    dgq(5, 13:15) =       q2_t(10:12)
    dgq(5, 16:18) =       q2_t(10:12)*this%phi(1, 2)
    dgq(5, 19:21) =       q2_t(10:12)*this%phi(2, 2)
    dgq(5, 22:24) =  ns + q2_t(10:12)*this%phi(3, 2)
    
    return
  end subroutine dgqrevolutejoint12

  subroutine krevolutejoint12(this, lambda, k, qopt)
    implicit none
    class(constraint_12_revolutejoint), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)

    k(:, :)      = 0.0d0
    
    k(1:3,1:3)   =   0.0d0
    k(1:3,4:6)   = - i_33*lambda(1)
    k(1:3,7:9)   = - i_33*lambda(2)
    k(1:3,10:12) = - i_33*lambda(3)
    k(1:3,16:18) = - i_33*lambda(4)
    k(1:3,22:24) = - i_33*lambda(5)
    
    k(4:6,1:3)   = - i_33*lambda(1)
    k(4:6,4:6)   =   i_33*(- this%phi(1,1)*lambda(1) - this%phi(1,1)*lambda(1))
    k(4:6,7:9)   =   i_33*(- this%phi(1,1)*lambda(2) - this%phi(2,1)*lambda(1))
    k(4:6,10:12) =   i_33*(- this%phi(1,1)*lambda(3) - this%phi(3,1)*lambda(1))
    k(4:6,13:15) =   i_33*lambda(1)
    k(4:6,16:18) =   i_33*(- this%phi(1,1)*lambda(4) + this%phi(1,2)*lambda(1))
    k(4:6,19:21) =   i_33*(                          + this%phi(2,2)*lambda(1))
    k(4:6,22:24) =   i_33*(- this%phi(1,1)*lambda(5) + this%phi(3,2)*lambda(1))
    
    k(7:9,1:3)   = - i_33*lambda(2)
    k(7:9,4:6)   = i_33*(- this%phi(2,1)*lambda(1) - this%phi(1,1)*lambda(2))
    k(7:9,7:9)   = i_33*(- this%phi(2,1)*lambda(2) - this%phi(2,1)*lambda(2))
    k(7:9,10:12) = i_33*(- this%phi(2,1)*lambda(3) - this%phi(3,1)*lambda(2))
    k(7:9,13:15) = i_33*lambda(2)
    k(7:9,16:18) = i_33*(- this%phi(2,1)*lambda(4) + this%phi(1,2)*lambda(2))
    k(7:9,19:21) = i_33*(                          + this%phi(2,2)*lambda(2))
    k(7:9,22:24) = i_33*(- this%phi(2,1)*lambda(5) + this%phi(3,2)*lambda(2))
    
    k(10:12,1:3)   = - i_33*lambda(3)
    k(10:12,4:6)   = i_33*(- this%phi(3,1)*lambda(1) - this%phi(1,1)*lambda(3))
    k(10:12,7:9)   = i_33*(- this%phi(3,1)*lambda(2) - this%phi(2,1)*lambda(3))
    k(10:12,10:12) = i_33*(- this%phi(3,1)*lambda(3) - this%phi(3,1)*lambda(3))
    k(10:12,13:15) = i_33*lambda(3)
    k(10:12,16:18) = i_33*(- this%phi(3,1)*lambda(4) + this%phi(1,2)*lambda(3))
    k(10:12,19:21) = i_33*(                          + this%phi(2,2)*lambda(3))
    k(10:12,22:24) = i_33*(- this%phi(3,1)*lambda(5) + this%phi(3,2)*lambda(3))
    
    k(13:15,4:6)   = i_33*lambda(1)
    k(13:15,7:9)   = i_33*lambda(2)
    k(13:15,10:12) = i_33*lambda(3)
    k(13:15,16:18) = i_33*lambda(4)
    k(13:15,22:24) = i_33*lambda(5)
    
    k(16:18,1:3)   = - i_33*lambda(4)
    k(16:18,4:6)   = i_33*(this%phi(1,2)*lambda(1) - this%phi(1,1)*lambda(4))
    k(16:18,7:9)   = i_33*(this%phi(1,2)*lambda(2) - this%phi(2,1)*lambda(4))
    k(16:18,10:12) = i_33*(this%phi(1,2)*lambda(3) - this%phi(3,1)*lambda(4))
    k(16:18,13:15) = i_33*lambda(4)
    k(16:18,16:18) = i_33*(this%phi(1,2)*lambda(4) + this%phi(1,2)*lambda(4))
    k(16:18,19:21) = i_33*(                        + this%phi(2,2)*lambda(4))
    k(16:18,22:24) = i_33*(this%phi(1,2)*lambda(5) + this%phi(3,2)*lambda(4))
    
    k(19:21,4:6)   = this%phi(2,2)*lambda(1)*i_33
    k(19:21,7:9)   = this%phi(2,2)*lambda(2)*i_33 
    k(19:21,10:12) = this%phi(2,2)*lambda(3)*i_33
    k(19:21,16:18) = this%phi(2,2)*lambda(4)*i_33
    k(19:21,22:24) = this%phi(2,2)*lambda(5)*i_33
    
    k(22:24,1:3)   = - i_33*lambda(5)
    k(22:24,4:6)   = i_33*(this%phi(3,2)*lambda(1) - this%phi(1,1)*lambda(5))
    k(22:24,7:9)   = i_33*(this%phi(3,2)*lambda(2) - this%phi(2,1)*lambda(5))
    k(22:24,10:12) = i_33*(this%phi(3,2)*lambda(3) - this%phi(3,1)*lambda(5))
    k(22:24,13:15) = i_33*lambda(5)
    k(22:24,16:18) = i_33*(this%phi(3,2)*lambda(4) + this%phi(1,2)*lambda(5))
    k(22:24,19:21) = i_33*(                        + this%phi(2,2)*lambda(5))
    k(22:24,22:24) = i_33*(this%phi(3,2)*lambda(5) + this%phi(3,2)*lambda(5))
    
    return
    end subroutine krevolutejoint12

! ====================================================================================================================
! ROTOR REVOLUTE JOINT
! ====================================================================================================================
! MIN 04/2024: A revolute joint constraint for rotating node B around axis prescribed by offset axes c and d to node A.
  subroutine grotorrevolutejoint12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_rotorrevolutejoint), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) :: d, e, d0, e0
    
    ! axes to which the angles of node B shall be constrained at t and t0
    d  = q2_t(1:3)-(q1_t(1:3)+this%phi(1,1)*q1_t(4:6)+this%phi(2,1)*q1_t(7:9)+this%phi(3,1)*q1_t(10:12))
    e  = q2_t(1:3)-(q1_t(1:3)+this%phi(1,2)*q1_t(4:6)+this%phi(2,2)*q1_t(7:9)+this%phi(3,2)*q1_t(10:12))
    
    d0 = q2_0(1:3)-(q1_0(1:3)+this%phi(1,1)*q1_0(4:6)+this%phi(2,1)*q1_0(7:9)+this%phi(3,1)*q1_0(10:12))
    e0 = q2_0(1:3)-(q1_0(1:3)+this%phi(1,2)*q1_0(4:6)+this%phi(2,2)*q1_0(7:9)+this%phi(3,2)*q1_0(10:12))
    
    ! constraints to maintain angles between axis d and directors at node B, and length of d
    g(1) = dot_product(q2_t( 4: 6),d) - dot_product(q2_0( 4: 6),d0)
    g(2) = dot_product(q2_t( 7: 9),d) - dot_product(q2_0( 7: 9),d0)
    g(3) = dot_product(q2_t(10:12),d) - dot_product(q2_0(10:12),d0)

    ! constraints to maintain angles between axis e and directors at node B, and length of e
    g(4) = dot_product(q2_t( 4: 6),e) - dot_product(q2_0( 4: 6),e0)
    g(5) = dot_product(q2_t( 7: 9),e) - dot_product(q2_0( 7: 9),e0)
    g(6) = dot_product(q2_t(10:12),e) - dot_product(q2_0(10:12),e0)

  end subroutine grotorrevolutejoint12

  subroutine dgqrotorrevolutejoint12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_rotorrevolutejoint), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(3) :: d, e
    
    d  = q2_t(1:3)-(q1_t(1:3)+this%phi(1,1)*q1_t(4:6)+this%phi(2,1)*q1_t(7:9)+this%phi(3,1)*q1_t(10:12))
    e  = q2_t(1:3)-(q1_t(1:3)+this%phi(1,2)*q1_t(4:6)+this%phi(2,2)*q1_t(7:9)+this%phi(3,2)*q1_t(10:12))
            
    dgq(:, :)     = 0.0d0
    
    dgq(1, 1: 3) = -q2_t(4:6)
    dgq(1, 4: 6) = -this%phi(1,1)*q2_t(4:6)
    dgq(1, 7: 9) = -this%phi(2,1)*q2_t(4:6)
    dgq(1,10:12) = -this%phi(3,1)*q2_t(4:6)
    dgq(1,13:15) =  q2_t(4:6)
    dgq(1,16:18) =  d
    
    dgq(2, 1: 3) = -q2_t(7:9)
    dgq(2, 4: 6) = -this%phi(1,1)*q2_t(7:9)
    dgq(2, 7: 9) = -this%phi(2,1)*q2_t(7:9)
    dgq(2,10:12) = -this%phi(3,1)*q2_t(7:9)
    dgq(2,13:15) =  q2_t(7:9)
    dgq(2,19:21) =  d
    
    dgq(3, 1: 3) = -q2_t(10:12)
    dgq(3, 4: 6) = -this%phi(1,1)*q2_t(10:12)
    dgq(3, 7: 9) = -this%phi(2,1)*q2_t(10:12)
    dgq(3,10:12) = -this%phi(3,1)*q2_t(10:12)
    dgq(3,13:15) =  q2_t(10:12)
    dgq(3,22:24) =  d
    
    dgq(4, 1: 3) = -q2_t(4:6)
    dgq(4, 4: 6) = -this%phi(1,2)*q2_t(4:6)
    dgq(4, 7: 9) = -this%phi(2,2)*q2_t(4:6)
    dgq(4,10:12) = -this%phi(3,2)*q2_t(4:6)
    dgq(4,13:15) =  q2_t(4:6)
    dgq(4,16:18) =  e
    
    dgq(5, 1: 3) = -q2_t(7:9)
    dgq(5, 4: 6) = -this%phi(1,2)*q2_t(7:9)
    dgq(5, 7: 9) = -this%phi(2,2)*q2_t(7:9)
    dgq(5,10:12) = -this%phi(3,2)*q2_t(7:9)
    dgq(5,13:15) =  q2_t(7:9)
    dgq(5,19:21) =  e
    
    dgq(6, 1: 3) = -q2_t(10:12)
    dgq(6, 4: 6) = -this%phi(1,2)*q2_t(10:12)
    dgq(6, 7: 9) = -this%phi(2,2)*q2_t(10:12)
    dgq(6,10:12) = -this%phi(3,2)*q2_t(10:12)
    dgq(6,13:15) =  q2_t(10:12)
    dgq(6,22:24) =  e
    
    return
  end subroutine dgqrotorrevolutejoint12

  subroutine krotorrevolutejoint12(this, lambda, k, qopt)
    implicit none
    class(constraint_12_rotorrevolutejoint), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)

    k(:, :)      = 0.0d0
    
    k( 1: 3,16:18) = -i_33*(lambda(1)+lambda(4))
    k( 1: 3,19:21) = -i_33*(lambda(2)+lambda(5))
    k( 1: 3,22:24) = -i_33*(lambda(3)+lambda(6))
    
    k( 4: 6,16:18) = -i_33*(this%phi(1,1)*lambda(1)+this%phi(1,2)*lambda(4))
    k( 4: 6,19:21) = -i_33*(this%phi(1,1)*lambda(2)+this%phi(1,2)*lambda(5))
    k( 4: 6,22:24) = -i_33*(this%phi(1,1)*lambda(3)+this%phi(1,2)*lambda(6))
    
    k( 7: 9,16:18) = -i_33*(this%phi(2,1)*lambda(1)+this%phi(2,2)*lambda(4))
    k( 7: 9,19:21) = -i_33*(this%phi(2,1)*lambda(2)+this%phi(2,2)*lambda(5))
    k( 7: 9,22:24) = -i_33*(this%phi(2,1)*lambda(3)+this%phi(2,2)*lambda(6))
    
    k(10:12,16:18) = -i_33*(this%phi(3,1)*lambda(1)+this%phi(3,2)*lambda(4))
    k(10:12,19:21) = -i_33*(this%phi(3,1)*lambda(2)+this%phi(3,2)*lambda(5))
    k(10:12,22:24) = -i_33*(this%phi(3,1)*lambda(3)+this%phi(3,2)*lambda(6))
    
    k(13:15,16:18) =  i_33*(lambda(1)+lambda(4))
    k(13:15,19:21) =  i_33*(lambda(2)+lambda(5))
    k(13:15,22:24) =  i_33*(lambda(3)+lambda(6))
    
    k(16:18, 1: 3) = -i_33*(lambda(1)+lambda(4))
    k(16:18, 4: 6) = -i_33*(this%phi(1,1)*lambda(1)+this%phi(1,2)*lambda(4))
    k(16:18, 7: 9) = -i_33*(this%phi(2,1)*lambda(1)+this%phi(2,2)*lambda(4))
    k(16:18,10:12) = -i_33*(this%phi(3,1)*lambda(1)+this%phi(3,2)*lambda(4))
    k(16:18,13:15) =  i_33*(lambda(1)+lambda(4))
    
    k(19:21, 1: 3) = -i_33*(lambda(2)+lambda(5))
    k(19:21, 4: 6) = -i_33*(this%phi(1,1)*lambda(2)+this%phi(1,2)*lambda(5))
    k(19:21, 7: 9) = -i_33*(this%phi(2,1)*lambda(2)+this%phi(2,2)*lambda(5))
    k(19:21,10:12) = -i_33*(this%phi(3,1)*lambda(2)+this%phi(3,2)*lambda(5))
    k(19:21,13:15) =  i_33*(lambda(2)+lambda(5))
    
    k(22:24, 1: 3) = -i_33*(lambda(3)+lambda(6))
    k(22:24, 4: 6) = -i_33*(this%phi(1,1)*lambda(3)+this%phi(1,2)*lambda(6))
    k(22:24, 7: 9) = -i_33*(this%phi(2,1)*lambda(3)+this%phi(2,2)*lambda(6))
    k(22:24,10:12) = -i_33*(this%phi(3,1)*lambda(3)+this%phi(3,2)*lambda(6))
    k(22:24,13:15) =  i_33*(lambda(3)+lambda(6))
    
    return
    end subroutine krotorrevolutejoint12

! ====================================================================================================================
! NEW REVOLUTE JOINT FROM GERADIN+CARDONA
! ====================================================================================================================
  subroutine grevolutejoint12_3(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_revolutejoint_3), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) :: ns, ns0
    
    !< conservation of the distance 
    g(1:3) = q2_t(1:3)-q1_t(1:3) - (q2_t(1:3)-q1_t(1:3))

    !< conservation of the angles
    g(4)   =  dot_product(q1_t(4:6), q2_t(10:12)) - dot_product(q1_0(4:6), q2_0(10:12))
    g(5)   =  dot_product(q1_t(7:9), q2_t(10:12)) - dot_product(q1_0(7:9), q2_0(10:12))

    return
    end subroutine grevolutejoint12_3
    
  subroutine dgqrevolutejoint12_3(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_revolutejoint_3), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(3) :: ns
            
    dgq(:, :)     = 0.0d0
    
    dgq(1:3,1: 3) = -i_33
    dgq(1:3,13:15)=  i_33
    
    dgq(4,  4: 6) = q2_t(10:12)
    dgq(4, 22:24) = q1_t(4:6)
    
    dgq(5 , 1: 3) = q2_t(10:12)
    dgq(5 ,22:24) = q1_t(7:9)

    return
  end subroutine dgqrevolutejoint12_3
    
subroutine krevolutejoint12_3(this, lambda, k, qopt)
    implicit none
    class(constraint_12_revolutejoint_3), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)

    k(:, :)      = 0.0d0

    k(4:6,22:24)   = i_33*lambda(4)
    k(7:9,22:24)   = i_33*lambda(5)
    k(19:21,7:9)   = i_33*lambda(4)
    k(22:24,10:12) = i_33*lambda(5)
    
    return    
  end subroutine krevolutejoint12_3
    
    
    
! ====================================================================================================================
! NEW REVOLUTE JOINT
! ====================================================================================================================
  subroutine grevolutejoint12_2(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_revolutejoint_2), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) :: ns, ns0
    
    !< rotation axis at time t
    ns  =   (q2_t(1:3)+this%phi(1,2)*q2_t(4:6)+this%phi(2,2)*q2_t(7:9)+this%phi(3,2)*q2_t(10:12)) & 
          - (q1_t(1:3)+this%phi(1,1)*q1_t(4:6)+this%phi(2,1)*q1_t(7:9)+this%phi(3,1)*q1_t(10:12))

    !< rotation axis at time t=0
    ns0 =   (q2_0(1:3)+this%phi(1, 2)*q2_0(4:6)+this%phi(2, 2)*q2_0(7:9)+this%phi(3, 2)*q2_0(10:12)) & 
          - (q1_0(1:3)+this%phi(1, 1)*q1_0(4:6)+this%phi(2, 1)*q1_0(7:9)+this%phi(3, 1)*q1_0(10:12))
  
    this%dir_gl = ns/norm2(ns) ! in global KOS to store for later use, if necessary

    !< conservation of angle and distance between distance vector at t and t0
    g(1) = dot_product(ns, q1_t(4:6))   - dot_product(ns0, q1_0(4:6))
    g(2) = dot_product(ns, q1_t(7:9))   - dot_product(ns0, q1_0(7:9))
    g(3) = dot_product(ns, q1_t(10:12)) - dot_product(ns0, q1_0(10:12))
    g(4) = dot_product(ns, q2_t(4:6))   - dot_product(ns0, q2_0(4:6))
    g(5) = dot_product(ns, q2_t(7:9))   - dot_product(ns0, q2_0(7:9))
    g(6) = dot_product(ns, q2_t(10:12)) - dot_product(ns0, q2_0(10:12))
    
    return
  end subroutine grevolutejoint12_2

  subroutine dgqrevolutejoint12_2(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_revolutejoint_2), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(3) :: ns

    ns  = + (q2_t(1:3)+this%phi(1, 2)*q2_t(4:6)+this%phi(2, 2)*q2_t(7:9)+this%phi(3, 2)*q2_t(10:12)) & 
          - (q1_t(1:3)+this%phi(1, 1)*q1_t(4:6)+this%phi(2, 1)*q1_t(7:9)+this%phi(3, 1)*q1_t(10:12))
            
    dgq(:, :)     = 0.0d0
    
    dgq(1,  1: 3) =     - q1_t(4:6)
    dgq(1,  4: 6) =  ns - q1_t(4:6)*this%phi(1, 1)
    dgq(1,  7: 9) =     - q1_t(4:6)*this%phi(2, 1)
    dgq(1, 10:12) =     - q1_t(4:6)*this%phi(3, 1)
    dgq(1, 13:15) =       q1_t(4:6)
    dgq(1, 16:18) =       q1_t(4:6)*this%phi(1, 2)
    dgq(1, 19:21) =       q1_t(4:6)*this%phi(2, 2)
    dgq(1, 22:24) =       q1_t(4:6)*this%phi(3, 2)
    
    dgq(2,  1: 3) =     - q1_t(7:9)
    dgq(2,  4: 6) =     - q1_t(7:9)*this%phi(1, 1)
    dgq(2,  7: 9) =  ns - q1_t(7:9)*this%phi(2, 1)
    dgq(2, 10:12) =     - q1_t(7:9)*this%phi(3, 1)
    dgq(2, 13:15) =       q1_t(7:9)
    dgq(2, 16:18) =       q1_t(7:9)*this%phi(1, 2)
    dgq(2, 19:21) =       q1_t(7:9)*this%phi(2, 2) 
    dgq(2, 22:24) =       q1_t(7:9)*this%phi(3, 2)
    
    dgq(3,  1: 3) =     - q1_t(10:12)
    dgq(3,  4: 6) =     - q1_t(10:12)*this%phi(1, 1)
    dgq(3,  7: 9) =     - q1_t(10:12)*this%phi(2, 1)
    dgq(3, 10:12) =  ns - q1_t(10:12)*this%phi(3, 1) 
    dgq(3, 13:15) =       q1_t(10:12)
    dgq(3, 16:18) =       q1_t(10:12)*this%phi(1, 2)
    dgq(3, 19:21) =       q1_t(10:12)*this%phi(2, 2) 
    dgq(3, 22:24) =       q1_t(10:12)*this%phi(3, 2)
    
    dgq(4,  1: 3) =     - q2_t(4:6)
    dgq(4,  4: 6) =     - q2_t(4:6)*this%phi(1, 1)
    dgq(4,  7: 9) =     - q2_t(4:6)*this%phi(2, 1)
    dgq(4, 10:12) =     - q2_t(4:6)*this%phi(3, 1)
    dgq(4, 13:15) =       q2_t(4:6)
    dgq(4, 16:18) =  ns + q2_t(4:6)*this%phi(1, 2)
    dgq(4, 19:21) =       q2_t(4:6)*this%phi(2, 2)
    dgq(4, 22:24) =       q2_t(4:6)*this%phi(3, 2)
    
    dgq(5,  1: 3) =     - q2_t(7:9)
    dgq(5,  4: 6) =     - q2_t(7:9)*this%phi(1, 1)
    dgq(5,  7: 9) =     - q2_t(7:9)*this%phi(2, 1)
    dgq(5, 10:12) =     - q2_t(7:9)*this%phi(3, 1)
    dgq(5, 13:15) =       q2_t(7:9)
    dgq(5, 16:18) =       q2_t(7:9)*this%phi(1, 2)
    dgq(5, 19:21) =  ns + q2_t(7:9)*this%phi(2, 2)
    dgq(5, 22:24) =       q2_t(7:9)*this%phi(3, 2)
    
    dgq(6,  1: 3) =     - q2_t(10:12)
    dgq(6,  4: 6) =     - q2_t(10:12)*this%phi(1, 1)
    dgq(6,  7: 9) =     - q2_t(10:12)*this%phi(2, 1)
    dgq(6, 10:12) =     - q2_t(10:12)*this%phi(3, 1)
    dgq(6, 13:15) =       q2_t(10:12)
    dgq(6, 16:18) =       q2_t(10:12)*this%phi(1, 2)
    dgq(6, 19:21) =       q2_t(10:12)*this%phi(2, 2)
    dgq(6, 22:24) =  ns + q2_t(10:12)*this%phi(3, 2)

    return
  end subroutine dgqrevolutejoint12_2

  subroutine krevolutejoint12_2(this, lambda, k, qopt)
    implicit none
    class(constraint_12_revolutejoint_2), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)

    k(:, :)      = 0.0d0
    
    k(1:3,4:6)   = - i_33*lambda(1)
    k(1:3,7:9)   = - i_33*lambda(2)
    k(1:3,10:12) = - i_33*lambda(3)
    k(1:3,16:18) = - i_33*lambda(4)
    k(1:3,19:21) = - i_33*lambda(5)
    k(1:3,22:24) = - i_33*lambda(6)
    
    k(4:6,1:3)   = - i_33*lambda(1)
    k(4:6,4:6)   = i_33*(- this%phi(1,1)*lambda(1) - this%phi(1,1)*lambda(1))
    k(4:6,7:9)   = i_33*(- this%phi(1,1)*lambda(2) - this%phi(2,1)*lambda(1))
    k(4:6,10:12) = i_33*(- this%phi(1,1)*lambda(3) - this%phi(3,1)*lambda(1))
    k(4:6,13:15) = i_33*lambda(1)
    k(4:6,16:18) = i_33*(- this%phi(1,1)*lambda(4) + this%phi(1,2)*lambda(1))
    k(4:6,19:21) = i_33*(- this%phi(1,1)*lambda(5) + this%phi(2,2)*lambda(1))
    k(4:6,22:24) = i_33*(- this%phi(1,1)*lambda(6) + this%phi(3,2)*lambda(1))
    
    k(7:9,1:3)   = - i_33*lambda(2)
    k(7:9,4:6)   = i_33*(- this%phi(2,1)*lambda(1) - this%phi(1,1)*lambda(2))
    k(7:9,7:9)   = i_33*(- this%phi(2,1)*lambda(2) - this%phi(2,1)*lambda(2))
    k(7:9,10:12) = i_33*(- this%phi(2,1)*lambda(3) - this%phi(3,1)*lambda(2))
    k(7:9,13:15) = i_33*lambda(2)
    k(7:9,16:18) = i_33*(- this%phi(2,1)*lambda(4) + this%phi(1,2)*lambda(2))
    k(7:9,19:21) = i_33*(- this%phi(2,1)*lambda(5) + this%phi(2,2)*lambda(2))
    k(7:9,22:24) = i_33*(- this%phi(2,1)*lambda(6) + this%phi(3,2)*lambda(2))
    
    k(10:12,1:3)   = - i_33*lambda(3)
    k(10:12,4:6)   = i_33*(- this%phi(3,1)*lambda(1) - this%phi(1,1)*lambda(3))
    k(10:12,7:9)   = i_33*(- this%phi(3,1)*lambda(2) - this%phi(2,1)*lambda(3))
    k(10:12,10:12) = i_33*(- this%phi(3,1)*lambda(3) - this%phi(3,1)*lambda(3))
    k(10:12,13:15) = i_33*lambda(3)
    k(10:12,16:18) = i_33*(- this%phi(3,1)*lambda(4) + this%phi(1,2)*lambda(3))
    k(10:12,19:21) = i_33*(- this%phi(3,1)*lambda(5) + this%phi(2,2)*lambda(3))
    k(10:12,22:24) = i_33*(- this%phi(3,1)*lambda(6) + this%phi(3,2)*lambda(3))
    
    k(13:15,4:6)   = i_33*lambda(1)
    k(13:15,7:9)   = i_33*lambda(2)
    k(13:15,10:12) = i_33*lambda(3)
    k(13:15,16:18) = i_33*lambda(4)
    k(13:15,19:21) = i_33*lambda(5)
    k(13:15,22:24) = i_33*lambda(6)
    
    k(16:18,1:3)   = - i_33*lambda(4)
    k(16:18,4:6)   = i_33*(this%phi(1,2)*lambda(1) - this%phi(1,1)*lambda(4))
    k(16:18,7:9)   = i_33*(this%phi(1,2)*lambda(2) - this%phi(2,1)*lambda(4))
    k(16:18,10:12) = i_33*(this%phi(1,2)*lambda(3) - this%phi(3,1)*lambda(4))
    k(16:18,13:15) = i_33*lambda(4)
    k(16:18,16:18) = i_33*(this%phi(1,2)*lambda(4) + this%phi(1,2)*lambda(4))
    k(16:18,19:21) = i_33*(this%phi(1,2)*lambda(5) + this%phi(2,2)*lambda(4))
    k(16:18,22:24) = i_33*(this%phi(1,2)*lambda(6) + this%phi(3,2)*lambda(4))
    
    k(19:21,1:3)   = - i_33*lambda(5)
    k(19:21,4:6)   = i_33*(this%phi(2,2)*lambda(1) - this%phi(1,1)*lambda(5))
    k(19:21,7:9)   = i_33*(this%phi(2,2)*lambda(2) - this%phi(2,1)*lambda(5))
    k(19:21,10:12) = i_33*(this%phi(2,2)*lambda(3) - this%phi(3,1)*lambda(5))
    k(19:21,13:15) = i_33*lambda(5)
    k(19:21,16:18) = i_33*(this%phi(2,2)*lambda(4) + this%phi(1,2)*lambda(5))
    k(19:21,19:21) = i_33*(this%phi(2,2)*lambda(5) + this%phi(2,2)*lambda(5))
    k(19:21,22:24) = i_33*(this%phi(2,2)*lambda(6) + this%phi(3,2)*lambda(5))  
    
    k(22:24,1:3)   = - i_33*lambda(6)
    k(22:24,4:6)   = i_33*(this%phi(3,2)*lambda(1) - this%phi(1,1)*lambda(6))
    k(22:24,7:9)   = i_33*(this%phi(3,2)*lambda(2) - this%phi(2,1)*lambda(6))
    k(22:24,10:12) = i_33*(this%phi(3,2)*lambda(3) - this%phi(3,1)*lambda(6))
    k(22:24,13:15) = i_33*lambda(6)
    k(22:24,16:18) = i_33*(this%phi(3,2)*lambda(4) + this%phi(1,2)*lambda(6))
    k(22:24,19:21) = i_33*(this%phi(3,2)*lambda(5) + this%phi(2,2)*lambda(6))
    k(22:24,22:24) = i_33*(this%phi(3,2)*lambda(6) + this%phi(3,2)*lambda(6))
    
    return    
  end subroutine krevolutejoint12_2
  
! ====================================================================================================================
! RIGID CONNECTION
! ====================================================================================================================  
  subroutine grigidconnection12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_rigidconnection), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: ns(3), ns0(3)
    integer :: i,j

    !< rotation axis at time t
    ns  =   (q2_t(1:3)+this%phi(1, 2)*q2_t(4:6)+this%phi(2, 2)*q2_t(7:9)+this%phi(3, 2)*q2_t(10:12)) & 
          - (q1_t(1:3)+this%phi(1, 1)*q1_t(4:6)+this%phi(2, 1)*q1_t(7:9)+this%phi(3, 1)*q1_t(10:12))
    
    !< rotation axis at time t=0
    ns0 =   (q2_0(1:3)+this%phi(1, 2)*q2_0(4:6)+this%phi(2, 2)*q2_0(7:9)+this%phi(3, 2)*q2_0(10:12)) & 
          - (q1_0(1:3)+this%phi(1, 1)*q1_0(4:6)+this%phi(2, 1)*q1_0(7:9)+this%phi(3, 1)*q1_0(10:12))
    
    !< conservation of angle between rotation axis and component ns on director of node 1
    g(1) = dot_product(ns,q1_t(4:6))   - dot_product(ns0,q1_0(4:6))
    g(2) = dot_product(ns,q1_t(7:9))   - dot_product(ns0,q1_0(7:9))
    g(3) = dot_product(ns,q1_t(10:12)) - dot_product(ns0,q1_0(10:12))

    do i = 1, 3
      do j = 1, 3
        g(i*3+j) = dot_product(q1_t(i*3+1:i*3+3), q2_t(j*3+1:j*3+3)) - &
                   dot_product(q1_0(i*3+1:i*3+3), q2_0(j*3+1:j*3+3))
      end do
    end do
    
    return
  end subroutine grigidconnection12

  subroutine dgqrigidconnection12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_rigidconnection), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8) :: ns(3), ns0(3)
    integer :: i,j
    
    ns  = + (q2_t(1:3)+this%phi(1, 2)*q2_t(4:6)+this%phi(2, 2)*q2_t(7:9)+this%phi(3, 2)*q2_t(10:12)) & 
          - (q1_t(1:3)+this%phi(1, 1)*q1_t(4:6)+this%phi(2, 1)*q1_t(7:9)+this%phi(3, 1)*q1_t(10:12))
            
    dgq(:, :)     = 0.0d0
    
    dgq(1,  1: 3) =     - q1_t(4:6)
    dgq(1,  4: 6) =  ns - q1_t(4:6)*this%phi(1, 1)
    dgq(1,  7: 9) =     - q1_t(4:6)*this%phi(2, 1)
    dgq(1, 10:12) =     - q1_t(4:6)*this%phi(3, 1)
    dgq(1, 13:15) =       q1_t(4:6)
    dgq(1, 16:18) =       q1_t(4:6)*this%phi(1, 2)
    dgq(1, 19:21) =       q1_t(4:6)*this%phi(2, 2)
    dgq(1, 22:24) =       q1_t(4:6)*this%phi(3, 2)
    
    dgq(2,  1: 3) =     - q1_t(7:9)
    dgq(2,  4: 6) =     - q1_t(7:9)*this%phi(1, 1)
    dgq(2,  7: 9) =  ns - q1_t(7:9)*this%phi(2, 1)
    dgq(2, 10:12) =     - q1_t(7:9)*this%phi(3, 1)
    dgq(2, 13:15) =       q1_t(7:9)
    dgq(2, 16:18) =       q1_t(7:9)*this%phi(1, 2)
    dgq(2, 19:21) =       q1_t(7:9)*this%phi(2, 2) 
    dgq(2, 22:24) =       q1_t(7:9)*this%phi(3, 2)
    
    dgq(3,  1: 3) =     - q1_t(10:12)
    dgq(3,  4: 6) =     - q1_t(10:12)*this%phi(1, 1)
    dgq(3,  7: 9) =     - q1_t(10:12)*this%phi(2, 1)
    dgq(3, 10:12) =  ns - q1_t(10:12)*this%phi(3, 1) 
    dgq(3, 13:15) =       q1_t(10:12)
    dgq(3, 16:18) =       q1_t(10:12)*this%phi(1, 2)
    dgq(3, 19:21) =       q1_t(10:12)*this%phi(2, 2) 
    dgq(3, 22:24) =       q1_t(10:12)*this%phi(3, 2)
    
    do i = 1, 3
      do j = 1, 3
        dgq(i*3+j, i*3+1 : i*3+3) = q2_t(j*3+1:j*3+3)
        dgq(i*3+j, j*3+13:j*3+15) = q1_t(i*3+1:i*3+3)
      end do
    end do
    
    return
  end subroutine dgqrigidconnection12

  subroutine krigidconnection12(this, lambda, k, qopt)
    implicit none
    class(constraint_12_rigidconnection), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    integer :: i,j

    k(:, :) = 0.0d0
    
    k(1:3,4:6)   = - i_33*lambda(1)
    k(1:3,7:9)   = - i_33*lambda(2)
    k(1:3,10:12) = - i_33*lambda(3)
    
    k(4:6,1:3)   = - i_33*lambda(1)
    k(4:6,4:6)   = i_33*(- this%phi(1,1)*lambda(1) - this%phi(1,1)*lambda(1))
    k(4:6,7:9)   = i_33*(- this%phi(1,1)*lambda(2) - this%phi(2,1)*lambda(1))
    k(4:6,10:12) = i_33*(- this%phi(1,1)*lambda(3) - this%phi(3,1)*lambda(1))
    k(4:6,13:15) = i_33*lambda(1)
    k(4:6,16:18) = i_33*(this%phi(1,2)*lambda(1))
    k(4:6,19:21) = i_33*(this%phi(2,2)*lambda(1))
    k(4:6,22:24) = i_33*(this%phi(3,2)*lambda(1))
    
    k(7:9,1:3)   = - i_33*lambda(2)
    k(7:9,4:6)   = i_33*(- this%phi(2,1)*lambda(1) - this%phi(1,1)*lambda(2))
    k(7:9,7:9)   = i_33*(- this%phi(2,1)*lambda(2) - this%phi(2,1)*lambda(2))
    k(7:9,10:12) = i_33*(- this%phi(2,1)*lambda(3) - this%phi(3,1)*lambda(2))
    k(7:9,13:15) = i_33*lambda(2)
    k(7:9,16:18) = i_33*(this%phi(1,2)*lambda(2))
    k(7:9,19:21) = i_33*(this%phi(2,2)*lambda(2))
    k(7:9,22:24) = i_33*(this%phi(3,2)*lambda(2))
    
    k(10:12,1:3)   = - i_33*lambda(3)
    k(10:12,4:6)   = i_33*(- this%phi(3,1)*lambda(1) - this%phi(1,1)*lambda(3))
    k(10:12,7:9)   = i_33*(- this%phi(3,1)*lambda(2) - this%phi(2,1)*lambda(3))
    k(10:12,10:12) = i_33*(- this%phi(3,1)*lambda(3) - this%phi(3,1)*lambda(3))
    k(10:12,13:15) = i_33*lambda(3)
    k(10:12,16:18) = i_33*(this%phi(1,2)*lambda(3))
    k(10:12,19:21) = i_33*(this%phi(2,2)*lambda(3))
    k(10:12,22:24) = i_33*(this%phi(3,2)*lambda(3))
    
    k(13:15,4:6)   = i_33*lambda(1)
    k(13:15,7:9)   = i_33*lambda(2)
    k(13:15,10:12) = i_33*lambda(3)
   
    k(16:18,4:6)   = i_33*(this%phi(1,2)*lambda(1))
    k(16:18,7:9)   = i_33*(this%phi(1,2)*lambda(2))
    k(16:18,10:12) = i_33*(this%phi(1,2)*lambda(3))
    
    k(19:21,4:6)   = i_33*(this%phi(2,2)*lambda(1))
    k(19:21,7:9)   = i_33*(this%phi(2,2)*lambda(2))
    k(19:21,10:12) = i_33*(this%phi(2,2)*lambda(3))
    
    k(22:24,4:6)   = i_33*(this%phi(3,2)*lambda(1))
    k(22:24,7:9)   = i_33*(this%phi(3,2)*lambda(2))
    k(22:24,10:12) = i_33*(this%phi(3,2)*lambda(3))
    
    do i = 1, 3
      do j = 1, 3
        k(i*3+1  : i*3+3,  j*3+13 : j*3+15) = k(i*3+1  : i*3+3,  j*3+13 : j*3+15) + lambda(i*3+j) * i_33
        k(j*3+13 : j*3+15, i*3+1  : i*3+3)  = k(j*3+13 : j*3+15, i*3+1  : i*3+3)  + lambda(i*3+j) * i_33 
      end do
    end do
    
    return
  end subroutine krigidconnection12
  
! ==================================================================================================================== 
! MASSPOINT CONNECTION
! ====================================================================================================================   
  subroutine gmasspointconnection12(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_masspointconnection), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: g(:)
    
    g = q1_t-q2_t
    
  end subroutine gmasspointconnection12

  subroutine dgqmasspointconnection12(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_masspointconnection), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)

    dgq(:, :) = 0.0d0
    dgq( 1:12,  1:12) = eye(12)
    dgq( 1:12, 13:24) =-eye(12)
    
  end subroutine dgqmasspointconnection12

! ====================================================================================================================   
! Angular velocity for a node12 around a defined fixed axis, regarding to fixed cos
! ====================================================================================================================   
  subroutine gangularvelocity_global(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_angularvelocity_global), intent(inout) :: this
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::  d1, d2, d3, d1_dot, d2_dot, d3_dot
    real(kind = 8), dimension(3) :: w_tilde, w
    
    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    d1_dot = (qopt1(4:6)   - q1_0(4:6))/this%deltat
    d2_dot = (qopt1(7:9)   - q1_0(7:9))/this%deltat
    d3_dot = (qopt1(10:12) - q1_0(10:12))/this%deltat
    
    w_tilde = this%dir*this%amplitude_rotation_vel
    w       = 0.5d0*( cross(d1,d1_dot) + cross(d2,d2_dot) + cross(d3,d3_dot))
    
    g(1) = dot_product(this%dir, w-w_tilde)

  end subroutine gangularvelocity_global

  subroutine dgqangularvelocity_global(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_angularvelocity_global), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::  d1, d2, d3
    
    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    dgq(:,:)     = 0.0d0
    dgq(1,4:6)   = 0.5d0*matmul(this%dir, skew(d1))
    dgq(1,7:9)   = 0.5d0*matmul(this%dir, skew(d2))
    dgq(1,10:12) = 0.5d0*matmul(this%dir, skew(d3))

  end subroutine dgqangularvelocity_global

  subroutine dgvangularvelocity_global(this, q1_t, q2_t, q1_0, q2_0, dgv, qopt1, qopt2)
    implicit none
    class(constraint_12_angularvelocity_global), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgv(:,:)
    real(kind = 8), dimension(3) ::  d1, d2, d3, d1_dot, d2_dot, d3_dot

    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    d1_dot = (qopt1(4:6)   - q1_0(4:6))/this%deltat
    d2_dot = (qopt1(7:9)   - q1_0(7:9))/this%deltat
    d3_dot = (qopt1(10:12) - q1_0(10:12))/this%deltat
    
    dgv(:,:) = 0.0d0
    dgv(1,4:6)   = 0.5d0*matmul(this%dir,skew(d1)/this%deltat - skew(d1_dot))
    dgv(1,7:9)   = 0.5d0*matmul(this%dir,skew(d2)/this%deltat - skew(d2_dot))
    dgv(1,10:12) = 0.5d0*matmul(this%dir,skew(d3)/this%deltat - skew(d3_dot))

  end subroutine dgvangularvelocity_global

  subroutine kangularvelocity_global(this, lambda, k, qopt)
    implicit none
    class(constraint_12_angularvelocity_global), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)

    k(:,:) = 0.0d0
    k(4:6,4:6)     = 0.5d0*skew(this%dir)*lambda(1)
    k(7:9,7:9)     = 0.5d0*skew(this%dir)*lambda(1)
    k(10:12,10:12) = 0.5d0*skew(this%dir)*lambda(1)
    
  end subroutine kangularvelocity_global

! ====================================================================================================================   
! Angular velocity for a node12 around a defined fixed axis, regarding to local cos
! ====================================================================================================================   
! g(qtn05(indices12a_q), qtn05(indices12b_q), qtn(indices12a_q), qtn(indices12b_q), constraint_g, qtn1(indices12a_q))
  subroutine gangularvelocity_local(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_angularvelocity_local), intent(inout) :: this
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::  d1, d2, d3, d1_dot, d2_dot, d3_dot
    real(kind = 8), dimension(3) :: w_tilde, w, n
    
    !< tn+1/2
    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    !< tn+1/2
    n = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3

! CHECK THIS ONE TOO!!!    
    d1_dot = (qopt1(4:6)   - q1_0(4:6))/this%deltat    ! d1_dot_{n+1/2} = (d1n+1-d1n)/delta
    d2_dot = (qopt1(7:9)   - q1_0(7:9))/this%deltat
    d3_dot = (qopt1(10:12) - q1_0(10:12))/this%deltat
    
    w_tilde = n*this%amplitude_rotation_vel
    w       = 0.5d0*( cross(d1,d1_dot) + cross(d2,d2_dot) + cross(d3,d3_dot))
    
    g(1) = dot_product(n, w-w_tilde)

  end subroutine gangularvelocity_local

  subroutine dgqangularvelocity_local(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_angularvelocity_local), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::  d1, d2, d3, n
    
    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)

    n = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3

    dgq(:,:)     = 0.0d0
    dgq(1,4:6)   = 0.5d0*matmul(n, skew(d1))
    dgq(1,7:9)   = 0.5d0*matmul(n, skew(d2))
    dgq(1,10:12) = 0.5d0*matmul(n, skew(d3))

  end subroutine dgqangularvelocity_local

  subroutine dgvangularvelocity_local(this, q1_t, q2_t, q1_0, q2_0, dgv, qopt1, qopt2)
    implicit none
    class(constraint_12_angularvelocity_local), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), allocatable, intent(inout) :: dgv(:,:)
    real(kind = 8), dimension(3) ::  d1, d2, d3, d1_dot, d2_dot, d3_dot, n
    real(kind = 8) :: A_n(3,9)
    
    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    d1_dot = (qopt1(4:6)   - q1_0(4:6))/this%deltat
    d2_dot = (qopt1(7:9)   - q1_0(7:9))/this%deltat
    d3_dot = (qopt1(10:12) - q1_0(10:12))/this%deltat

    n = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3
    
    A_n(:,:)     = 0.0d0
    A_n(1:3,1:3) = this%dir(1)*i_33
    A_n(1:3,4:6) = this%dir(2)*i_33
    A_n(1:3,7:9) = this%dir(3)*i_33
    
    dgv(:,:)     = 0.0d0
    dgv(1,4:6)   = 0.5d0*matmul(n,skew(d1)/this%deltat - skew(d1_dot)) + 0.5d0*matmul(d1,matmul(skew(d1_dot),A_n(1:3,1:3)))
    dgv(1,7:9)   = 0.5d0*matmul(n,skew(d2)/this%deltat - skew(d2_dot)) + 0.5d0*matmul(d2,matmul(skew(d2_dot),A_n(1:3,4:6)))
    dgv(1,10:12) = 0.5d0*matmul(n,skew(d3)/this%deltat - skew(d3_dot)) + 0.5d0*matmul(d3,matmul(skew(d3_dot),A_n(1:3,7:9)))

  end subroutine dgvangularvelocity_local

  subroutine kangularvelocity_local(this, lambda, k, qopt)
    implicit none
    class(constraint_12_angularvelocity_local), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    real(kind = 8), dimension(3) ::  d1, d2, d3, n
    real(kind = 8) :: A_n(3,9)
    
    d1 = qopt(4:6)
    d2 = qopt(7:9)
    d3 = qopt(10:12)
    
    n = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3
    
    A_n(:,:)     = 0.0d0
    A_n(1:3,1:3) = this%dir(1)*i_33
    A_n(1:3,4:6) = this%dir(2)*i_33
    A_n(1:3,7:9) = this%dir(3)*i_33
        
    k(:,:)         = 0.0d0
    k(4:6,4:6)     = 0.5d0*(skew(n) - matmul(skew(d1),A_n(1:3,1:3)))*lambda(1)
    k(7:9,7:9)     = 0.5d0*(skew(n) - matmul(skew(d2),A_n(1:3,4:6)))*lambda(1)
    k(10:12,10:12) = 0.5d0*(skew(n) - matmul(skew(d3),A_n(1:3,7:9)))*lambda(1)
    
    end subroutine kangularvelocity_local

! ====================================================================================================================   
! fixed or prescribed rotation for a node12 around a defined (corotational) axis, relative to a secend node, related to local cos
! ====================================================================================================================      
  subroutine grelative_rotation_local(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_relative_rotation_local), intent(inout) :: this
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::delx, delx_0
    real(kind = 8) :: phi
    
    delx   = q2_t(1:3)-q1_t(1:3)
    delx_0 = q2_0(1:3)-q1_0(1:3)
    
    phi = this%dir(1)
    
    g(1)     = dot_product(delx,q1_t(4:6))-dot_product(delx_0,q1_0(4:6))
    g(2)     = dot_product(delx,q1_t(7:9))-dot_product(delx_0,q1_0(7:9))
    g(3)     = dot_product(delx,q1_t(10:12))-dot_product(delx_0,q1_0(10:12))
    g(4:6)   = q2_t(4:6)-cos(phi)*q1_t(4:6)-sin(phi)*q1_t(7:9)
    g(7:9)   = q2_t(7:9)-cos(phi)*q1_t(7:9)-sin(phi)*q1_t(4:6)
    g(10:12) = q2_t(10:12)-q1_t(10:12)
    
    
  end subroutine grelative_rotation_local

  subroutine dgqrelative_rotation_local(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_relative_rotation_local), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    integer :: i,j,k,l
    real(kind = 8) :: phi
    
    phi = this%dir(1)
    
    dgq(:,:) = 0.0d0
    
    dgq(1,1:3) = -q1_t(4:6)
    dgq(1,4:6) = q2_t(1:3)-q1_t(1:3)
    dgq(1,13:15) = q1_t(4:6)
    
    dgq(2,1:3) = -q1_t(7:9)
    dgq(2,7:9) = q2_t(1:3)-q1_t(1:3)
    dgq(2,13:15) = q1_t(7:9)
    
    dgq(3,1:3) = -q1_t(10:12)
    dgq(3,10:12) = q2_t(1:3)-q1_t(1:3)
    dgq(3,13:15) = q1_t(10:12)
    
    dgq(4:6,4:6) = -cos(phi)*i_33
    dgq(4:6,7:9) = -sin(phi)*i_33
    dgq(4:6,16:18) = i_33
    
    dgq(7:9,4:6) = sin(phi)*i_33
    dgq(7:9,7:9) = -cos(phi)*i_33
    dgq(7:9,19:21) = i_33
    
    dgq(10:12,10:12) = -i_33
    dgq(10:12,22:24) = i_33           
    
    end subroutine dgqrelative_rotation_local

  subroutine krelative_rotation_local(this, lambda, k, qopt)
    implicit none
    class(constraint_12_relative_rotation_local), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    
    k(:,:)      = 0.0d0
    
    k(1:3,4:6)     = -lambda(1)*i_33
    k(1:3,7:9)     = -lambda(2)*i_33
    k(1:3,10:12)   = -lambda(3)*i_33
 
    k(4:6,1:3)     = -lambda(1)*i_33
    k(4:6,13:15)   =  lambda(1)*i_33
    
    k(7:9,1:3)     = -lambda(2)*i_33
    k(7:9,13:15)   =  lambda(2)*i_33
    
    k(10:12,1:3)   = -lambda(3)*i_33
    k(10:12,13:15) =  lambda(3)*i_33
    
    k(13:15,4:6)   =  lambda(1)*i_33
    k(13:15,7:9)   =  lambda(2)*i_33
    k(13:15,10:12) =  lambda(3)*i_33
        
 end subroutine krelative_rotation_local
    
! ====================================================================================================================   
! fixed or prescribed rotation for a node12 around a defined (corotational) axis, related to local cos
! ====================================================================================================================   
  subroutine grotation_local(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_rotation_local), intent(inout) :: this
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::  d1, d2, d3, delta_d1, delta_d2, delta_d3
    real(kind = 8), dimension(3) :: delta_phi, delta_phi_tilde, n
    
    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    delta_d1 = q1_t(4:6)   - qopt1(4:6)
    delta_d2 = q1_t(7:9)   - qopt1(7:9)
    delta_d3 = q1_t(10:12) - qopt1(10:12)

    delta_phi(:)       = 0.5d0*cross(d1,delta_d1) + 0.5d0*cross(d2,delta_d2) + 0.5d0*cross(d3,delta_d3)
    n = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3
    delta_phi_tilde(:) = this%amplitude_rotation*n
    
    g(1) = dot_product(n,delta_phi - delta_phi_tilde)

  end subroutine grotation_local

  subroutine dgqrotation_local(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_rotation_local), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::  d1, d2, d3, delta_d1, delta_d2, delta_d3, delta_phi, delta_phi_tilde, n
    real(kind = 8), dimension(1,12) :: A_dn(3,9), A_phi_tilde(3,9), A_n(3,9)

    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    delta_d1 = q1_t(4:6)   - qopt1(4:6)
    delta_d2 = q1_t(7:9)   - qopt1(7:9)
    delta_d3 = q1_t(10:12) - qopt1(10:12)

    delta_phi(:)       = 0.5d0*cross(d1,delta_d1) + 0.5d0*cross(d2,delta_d2) + 0.5d0*cross(d3,delta_d3)
    n(:)               = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3
    delta_phi_tilde(:) = this%amplitude_rotation*n
    
    A_n(:,:)     = 0.0d0
    A_n(1:3,1:3) = this%dir(1)*i_33
    A_n(1:3,4:6) = this%dir(2)*i_33
    A_n(1:3,7:9) = this%dir(3)*i_33
    
    A_dn(1:3,1:3) = 0.5d0*skew(qopt1(4:6))
    A_dn(1:3,4:6) = 0.5d0*skew(qopt1(7:9))
    A_dn(1:3,7:9) = 0.5d0*skew(qopt1(10:12))
    
    dgq(:,:)    = 0.0d0
    dgq(1,4:12) = matmul(delta_phi - delta_phi_tilde,A_n) + matmul(n,A_dn)    

  end subroutine dgqrotation_local

  subroutine krotation_local(this, lambda, k, qopt)
    implicit none
    class(constraint_12_rotation_local), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt(:)
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    real(kind = 8), dimension(1,12) :: A_dn(3,9), A_phi_tilde(3,9), A_n(3,9)

    A_n(:,:)     = 0.0d0
    A_n(1:3,1:3) = this%dir(1)*i_33
    A_n(1:3,4:6) = this%dir(2)*i_33
    A_n(1:3,7:9) = this%dir(3)*i_33
    
    A_dn(:,:)     = 0.0d0
    A_dn(1:3,1:3) = 0.5d0*skew(qopt(4:6))
    A_dn(1:3,4:6) = 0.5d0*skew(qopt(7:9))
    A_dn(1:3,7:9) = 0.5d0*skew(qopt(10:12))
        
    k(:,:)       = 0.0d0
    k(4:12,4:12) = lambda(1)*( matmul(transpose(A_dn), A_n) + matmul(transpose(A_n), A_dn) )
    
  end subroutine krotation_local

! ====================================================================================================================   
! fixed or prescribed rotation for a node12 around a defined axis, related to global cos
! ====================================================================================================================   
  subroutine grotation_global(this, q1_t, q2_t, q1_0, q2_0, g, qopt1, qopt2)
    implicit none
    class(constraint_12_rotation_global), intent(inout) :: this
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8), dimension(3) ::  d1, d2, d3, delta_d1, delta_d2, delta_d3
    real(kind = 8), dimension(3) :: delta_phi
    
    d1 = q1_t(4:6)
    d2 = q1_t(7:9)
    d3 = q1_t(10:12)
    
    delta_d1 = q1_t(4:6)   - qopt1(4:6)
    delta_d2 = q1_t(7:9)   - qopt1(7:9)
    delta_d3 = q1_t(10:12) - qopt1(10:12)

    delta_phi     = 0.5d0*( cross(d1,delta_d1) + cross(d2,delta_d2) + cross(d3,delta_d3)  )
    
    g(1) = dot_product(this%dir,delta_phi - this%amplitude_rotation*this%dir)

  end subroutine grotation_global

  subroutine dgqrotation_global(this, q1_t, q2_t, q1_0, q2_0, dgq, qopt1, qopt2)
    implicit none
    class(constraint_12_rotation_global), intent(in) :: this
    real(kind = 8), dimension(12), intent(in) :: q1_t, q2_t, q1_0, q2_0
    real(kind = 8), allocatable, intent(inout) :: dgq(:,:)
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: A_dn(3,9)
    
    A_dn(:,:)     = 0.0d0
    A_dn(1:3,1:3) = 0.5d0*skew(qopt1(4:6))
    A_dn(1:3,4:6) = 0.5d0*skew(qopt1(7:9))
    A_dn(1:3,7:9) = 0.5d0*skew(qopt1(10:12))

    dgq(:,:)    = 0.0d0
    dgq(1,4:12) = matmul(this%dir,A_dn)
    
  end subroutine dgqrotation_global
  
end module class_constraint_12
