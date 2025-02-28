module class_constraint_6

  use class_constraint
  use my_math_structure, only : outer, skew, eye
  use my_constants_structure, only : eps
  implicit none

  !> shell constraint base class
  type, extends(constraint) :: constraint_6   
  contains
     
     procedure :: gconstraint6 => gdummy6
     procedure :: dgconstraint6 => dgdummy6
     procedure :: kconstraint6 => kdummy6
  
  end type constraint_6
  
  type constraint_6_c   !! container type                                               
    class(constraint_6), pointer :: c
  end type
  
  !> shell internal constraint class
  type, extends(constraint_6) :: constraint_6_internal
  contains
     
    procedure :: gconstraint6 => ginternal6
    procedure :: dgconstraint6 => dginternal6
    procedure :: kconstraint6 => kinternal6
    
  end type constraint_6_internal
  
  interface constraint_6_internal
    module procedure create_constraint_6_internal
  end interface

!> shell rotation global constraint class
  type, extends(constraint_6) :: constraint_6_rotation_global
  contains
     
    procedure :: gconstraint6 => grotation_global6
    procedure :: dgconstraint6 => dgrotation_global6
    
  end type constraint_6_rotation_global
  
  interface constraint_6_rotation_global
    module procedure create_constraint_6_rotation_global
  end interface
  
  !> shell spherical support constraint class
  type, extends(constraint_6) :: constraint_6_sphericalsupport
  contains
     
    procedure :: gconstraint6 => gsphericalsupport6
    procedure :: dgconstraint6 => dgsphericalsupport6
    
  end type constraint_6_sphericalsupport
  
  interface constraint_6_sphericalsupport
    module procedure create_constraint_6_sphericalsupport
  end interface
     
  !> shell revolute support constraint class
  type, extends(constraint_6) :: constraint_6_revolutesupport
  contains
     
    procedure :: gconstraint6 => grevolutesupport6
    procedure :: dgconstraint6 => dgrevolutesupport6
    
  end type constraint_6_revolutesupport
  
  interface constraint_6_revolutesupport
    module procedure create_constraint_6_revolutesupport
  end interface

  !> shell spherical joint constraint class
  type, extends(constraint_6) :: constraint_6_sphericaljoint
  contains
     
    procedure :: gconstraint6 => gsphericaljoint6
    procedure :: dgconstraint6 => dgsphericaljoint6
    
  end type constraint_6_sphericaljoint
  
  interface constraint_6_sphericaljoint
    module procedure create_constraint_6_sphericaljoint
  end interface
  
  !> shell revolute joint constraint class
  type, extends(constraint_6) :: constraint_6_revolutejoint
  contains
     
    procedure :: gconstraint6 => grevolutejoint6
    procedure :: dgconstraint6 => dgrevolutejoint6
    
  end type constraint_6_revolutejoint
  
  interface constraint_6_revolutejoint
    module procedure create_constraint_6_revolutejoint
  end interface
  
  !> shell inextensible revolute joint constraint class
  type, extends(constraint_6) :: constraint_6_inextensiblerevolutejoint
  contains
     
    procedure :: gconstraint6 => ginextensiblerevolutejoint6
    procedure :: dgconstraint6 => dginextensiblerevolutejoint6
    procedure :: kconstraint6 => kinextensiblerevolutejoint6
    
  end type constraint_6_inextensiblerevolutejoint
  
  interface constraint_6_inextensiblerevolutejoint
    module procedure create_constraint_6_inextensiblerevolutejoint
  end interface
  
  !> shell layer connection constraint class
  type, extends(constraint_6) :: constraint_6_layerconnection
  contains
     
    procedure :: gconstraint6 => glayerconnection6
    procedure :: dgconstraint6 => dglayerconnection6
    
  end type constraint_6_layerconnection
  
  interface constraint_6_layerconnection
    module procedure create_constraint_6_layerconnection
  end interface
  
  !> shell simple support constraint class
  type, extends(constraint_6) :: constraint_6_simplesupport
  contains
     
    procedure :: gconstraint6 => gsimplesupport6
    procedure :: dgconstraint6 => dgsimplesupport6
    
  end type constraint_6_simplesupport
  
  interface constraint_6_simplesupport
    module procedure create_constraint_6_simplesupport
  end interface
  
  !> shell simple support constraint class
  type, extends(constraint_6) :: constraint_6_symmetry
  contains
     
    procedure :: gconstraint6 => gsymmetry6
    procedure :: dgconstraint6 => dgsymmetry6
    
  end type constraint_6_symmetry
  
  interface constraint_6_symmetry
    module procedure create_constraint_6_symmetry
  end interface

  contains
  
  subroutine gdummy6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    implicit none

    class(constraint_6), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: g(:)
    
  end subroutine gdummy6

  subroutine dgdummy6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none

    class(constraint_6), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    
  end subroutine dgdummy6

  subroutine kdummy6(this, lambda, k)

    implicit none

    class(constraint_6), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)

    real(kind = 8), allocatable, intent(inout) :: k(:,:)

  end subroutine kdummy6
  
  type(constraint_6_internal) function create_constraint_6_internal()
    
    create_constraint_6_internal%rankcount = 6
    create_constraint_6_internal%coordinates = 1
    create_constraint_6_internal%stiffness = .TRUE.
    return
  end function
  
  subroutine ginternal6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    
    implicit none

    class(constraint_6_internal), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: g(:)

    g = 0.5d0*(dot_product(q1_t(4:6), q1_t(4:6))-dot_product(q1_0(4:6), q1_0(4:6)))
    
  end subroutine ginternal6

  subroutine dginternal6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none

    class(constraint_6_internal), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)

    dg(1,:) = 0.0d0
    
    dg(1, 4:6) = q1_t(4:6)
    
  end subroutine dginternal6

  subroutine kinternal6(this, lambda, k)

    implicit none

    class(constraint_6_internal), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)

    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    
    k(:, :) = 0.0d0

    k(4:6, 4:6) = lambda(1)*i_33

  end subroutine kinternal6

  type(constraint_6_sphericalsupport) function create_constraint_6_sphericalsupport()
    
    create_constraint_6_sphericalsupport%rankcount = 3
    create_constraint_6_sphericalsupport%coordinates = 1
    create_constraint_6_sphericalsupport%stiffness = .FALSE.
    return
  end function
  
  subroutine gsphericalsupport6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)

    implicit none

    class(constraint_6_sphericalsupport), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: g(:)

    g = q1_t(1:3)-q1_0(1:3)

    return
    
  end subroutine gsphericalsupport6

  subroutine dgsphericalsupport6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)

    implicit none

    class(constraint_6_sphericalsupport), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)

    dg(:, :) = 0.0d0
    
    dg(1:3, 1:3) = i_33

    return
    
  end subroutine dgsphericalsupport6
  
  type(constraint_6_revolutesupport) function create_constraint_6_revolutesupport()
    
    create_constraint_6_revolutesupport%rankcount = 6
    create_constraint_6_revolutesupport%coordinates = 1
    create_constraint_6_revolutesupport%stiffness = .FALSE.
    return
  end function
  
  subroutine grevolutesupport6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    
    implicit none

    class(constraint_6_revolutesupport), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
    
    real(kind = 8), allocatable, intent(inout) :: g(:)
    
    g = q1_t-q1_0
    
    return
    
  end subroutine grevolutesupport6

  subroutine dgrevolutesupport6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none

    class(constraint_6_revolutesupport), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
        
    real(kind = 8), allocatable, intent(inout) :: dg(:,:)

    dg(:, :) = 0.0d0
    
    dg(1:3, 1:3) = i_33
    
    dg(4:6, 4:6) = i_33
    
    return
    
  end subroutine dgrevolutesupport6

  type(constraint_6_sphericaljoint) function create_constraint_6_sphericaljoint()
    
    create_constraint_6_sphericaljoint%rankcount = 3
    create_constraint_6_sphericaljoint%coordinates = 2
    create_constraint_6_sphericaljoint%stiffness = .FALSE.
    return
  end function
  
  subroutine gsphericaljoint6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    
    implicit none

    class(constraint_6_sphericaljoint), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
        
    real(kind = 8), allocatable, intent(inout) :: g(:)
        
    g(1:3) = q1_t(1:3)-q2_t(1:3)
    
    return
    
  end subroutine gsphericaljoint6

  subroutine dgsphericaljoint6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none

    class(constraint_6_sphericaljoint), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
        
    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    
    dg(:, :) = 0.0d0

    dg(1:3, 1:3) = i_33
    dg(1:3, 7:9) =-i_33
    
    return
    
  end subroutine dgsphericaljoint6

  type(constraint_6_revolutejoint) function create_constraint_6_revolutejoint()
    
    create_constraint_6_revolutejoint%rankcount = 6
    create_constraint_6_revolutejoint%coordinates = 2
    create_constraint_6_revolutejoint%stiffness = .FALSE.
    return
  end function
  
  subroutine grevolutejoint6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    
    implicit none

    class(constraint_6_revolutejoint), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: g(:)
    
    g(1:3) = q1_t(1:3)-q2_t(1:3)

    g(4:6) = q1_t(4:6)-q2_t(4:6)
    
    return
    
  end subroutine grevolutejoint6

  subroutine dgrevolutejoint6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none

    class(constraint_6_revolutejoint), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    
    dg(:, :) = 0.0d0

    dg(1:3, 1:3) = i_33
    dg(1:3, 7:9) =-i_33
    
    dg(4:6,  4: 6) = i_33
    dg(4:6, 10:12) =-i_33
    
    return
    
  end subroutine dgrevolutejoint6
  
  type(constraint_6_inextensiblerevolutejoint) function create_constraint_6_inextensiblerevolutejoint()
    create_constraint_6_inextensiblerevolutejoint%rankcount = 6
    create_constraint_6_inextensiblerevolutejoint%coordinates = 2
    create_constraint_6_inextensiblerevolutejoint%stiffness = .TRUE.
    return
  end function

  subroutine ginextensiblerevolutejoint6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    
    implicit none
  
    class(constraint_6_inextensiblerevolutejoint), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
  
    real(kind = 8), allocatable, intent(inout) :: g(:)
    
    g(1:3) = q1_t(1:3)-q2_t(1:3)
  
    g(4) = 0.5d0*(dot_product(q1_t(4:6), q1_t(4:6))-dot_product(q1_0(4:6), q1_0(4:6)))
    
    g(5) = 0.5d0*(dot_product(q2_t(4:6), q2_t(4:6))-dot_product(q2_0(4:6), q2_0(4:6)))
  
    g(6) = dot_product(q1_t(4:6), q2_t(4:6))-dot_product(q1_0(4:6), q2_0(4:6))
    
    return
    
  end subroutine ginextensiblerevolutejoint6
  
  subroutine dginextensiblerevolutejoint6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none
  
    class(constraint_6_inextensiblerevolutejoint), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
  
    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    
    dg(:, :) = 0.0d0
  
    ! spherical joint
    dg(1:3, 1:3) = i_33
    dg(1:3, 7:9) =-i_33
  
    ! internal
    dg(4, 4:6) = q1_t(4:6)
    dg(5, 10: 12) = q2_t(4:6)
    
    dg(6,  4: 6) = q2_t(4:6)
    dg(6, 10:12) = q1_t(4:6)
    
    return
    
  end subroutine dginextensiblerevolutejoint6
  
  subroutine kinextensiblerevolutejoint6(this, lambda, k)
    
    implicit none
    
    class(constraint_6_inextensiblerevolutejoint), intent(in) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)
    
    real(kind = 8), allocatable, intent(inout) :: k(:,:)
    
    k(:, :) = 0.0d0
    
    k( 4: 6,  4: 6) = lambda(4)*i_33
    k( 4: 6, 10:12) = lambda(6)*i_33
    k(10:12,  4: 6) = lambda(6)*i_33
    k(10:12, 10:12) = lambda(5)*i_33
    
    return
    
  end subroutine kinextensiblerevolutejoint6

  type(constraint_6_layerconnection) function create_constraint_6_layerconnection()
    
    create_constraint_6_layerconnection%rankcount = 3
    create_constraint_6_layerconnection%coordinates = 2
    create_constraint_6_layerconnection%stiffness = .FALSE.
    return
  end function
  
  subroutine glayerconnection6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    
    implicit none

    class(constraint_6_layerconnection), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: g(:)
    
    g(1:3) = q1_t(1:3)+0.5d0*this%phi(3, 1)*q1_t(4:6)-q2_t(1:3)+0.5d0*this%phi(3, 2)*q2_t(4:6)
    
    return
    
  end subroutine glayerconnection6

  subroutine dglayerconnection6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none

    class(constraint_6_layerconnection), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    
    dg(:, :) = 0.0d0

    dg(1:3,  1: 3) =                      i_33
    dg(1:3,  4: 6) = 0.5d0*this%phi(3, 1)*i_33
    dg(1:3,  7: 9) =-                     i_33
    dg(1:3, 10:12) = 0.5d0*this%phi(3, 2)*i_33   
    
    return
    
  end subroutine dglayerconnection6

  type(constraint_6_simplesupport) function create_constraint_6_simplesupport()
    
    create_constraint_6_simplesupport%rankcount = 1
    create_constraint_6_simplesupport%coordinates = 1
    create_constraint_6_simplesupport%stiffness = .FALSE.
    return
  end function
  
  subroutine gsimplesupport6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    
    implicit none

    class(constraint_6_simplesupport), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8), dimension(3) :: u_n, u_tilde
    
    u_n = q1_t(1:3) - q1_0(1:3)
    u_tilde = this%amplitude_translation*this%dir
    
    g = dot_product(u_n - u_tilde, this%dir)
    return
    
  end subroutine gsimplesupport6

  subroutine dgsimplesupport6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    
    implicit none

    class(constraint_6_simplesupport), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    
    dg(1,:) = 0.0d0
    dg(1, 1: 3) = this%dir

    return
    
  end subroutine dgsimplesupport6
  
  type(constraint_6_symmetry) function create_constraint_6_symmetry()
    
    create_constraint_6_symmetry%rankcount = 2
    create_constraint_6_symmetry%coordinates = 1
    create_constraint_6_symmetry%stiffness = .FALSE.
    return
  end function
  
  subroutine gsymmetry6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)

    implicit none

    class(constraint_6_symmetry), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
    
    real(kind = 8), allocatable, intent(inout) :: g(:)

    g(1) = dot_product(this%dir, q1_t(1:3)-q1_0(1:3))

    g(2) = dot_product(this%dir, q1_t(4:6))
    
    return
    
  end subroutine gsymmetry6

  subroutine dgsymmetry6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)

    implicit none

    class(constraint_6_symmetry), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)

    dg(:, :) = 0.0d0

    dg(1, 1:3) = this%dir

    dg(2, 4:6) = this%dir
    
    return
    
  end subroutine dgsymmetry6
  
  type(constraint_6_rotation_global) function create_constraint_6_rotation_global()
    create_constraint_6_rotation_global%rankcount = 3
    create_constraint_6_rotation_global%coordinates = 1
    create_constraint_6_rotation_global%stiffness = .FALSE.
    return
  end function
  
  subroutine grotation_global6(this, q1_t, q2_t, q1_0, q2_0, g, qopt1)
    implicit none
    class(constraint_6_rotation_global), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
    real(kind = 8), allocatable, intent(inout) :: g(:)
    real(kind = 8) :: R(3,3), amplitude_rotation
    
    amplitude_rotation = this%amplitude_rotation
    if (norm2(this%dir) .lt. eps) amplitude_rotation = 0.0d0
    R = dcos(amplitude_rotation) * i_33 + dsin(amplitude_rotation) * skew(this%dir) + (1.0d0-dcos(amplitude_rotation)) * outer(this%dir, this%dir)
    
    g(1:3) = q1_t(4:6) - matmul(R,q1_0(4:6))
    
  end subroutine grotation_global6

  subroutine dgrotation_global6(this, q1_t, q2_t, q1_0, q2_0, dg, qopt1)
    implicit none
    class(constraint_6_rotation_global), intent(in) :: this
    real(kind = 8), dimension(6), intent(in) :: q1_t, q2_t, q1_0, q2_0, qopt1
    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    real(kind = 8) :: d(3)
    
    dg(:,:)      = 0.0d0
    dg(1:3, 4:6) = eye(3)
    
  end subroutine dgrotation_global6
  
end module class_constraint_6
  
