module class_constraint_6_to_12

  use class_constraint

  implicit none

  type, extends(constraint) :: constraint_6to12

    contains

      procedure :: gconstraint6to12  => gdummy6to12
      procedure :: dgconstraint6to12 => dgdummy6to12
      procedure :: kconstraint6to12  => kdummy6to12
     
  end type constraint_6to12
   
  type constraint_6to12_c   !! container type                                               
    class(constraint_6to12), pointer :: c
  end type
  
  !> beam shell connector rigid cross section constraint
  type, extends(constraint_6to12) :: constraint_6to12_rigidtransition
  contains
     
    procedure :: gconstraint6to12  => grigidtransition6to12
    procedure :: dgconstraint6to12 => dgrigidtransition6to12
    procedure :: kconstraint6to12  => krigidtransition6to12
    
  end type constraint_6to12_rigidtransition
  
  interface constraint_6to12_rigidtransition
    module procedure create_constraint_6to12_rigidtransition
  end interface
     
  !> beam shell connector deformable cross section constraint
  type, extends(constraint_6to12) :: constraint_6to12_softtransition
  contains
     
    procedure :: gconstraint6to12  => gsofttransition6to12
    procedure :: dgconstraint6to12 => dgsofttransition6to12
    procedure :: kconstraint6to12  => ksofttransition6to12
    
  end type constraint_6to12_softtransition
  
  interface constraint_6to12_softtransition
    module procedure create_constraint_6to12_softtransition
  end interface
  
  
  contains
  
  !> q12 is the node12 side, q6 is the node6 side
  subroutine gdummy6to12(this, q12_t, q6_t, q12_0, q6_0, g)
    implicit none

    class(constraint_6to12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q12_t, q12_0
    real(kind = 8), dimension( 6), intent(in) :: q6_t, q6_0

    real(kind = 8), allocatable, intent(inout) :: g(:)
    
  end subroutine gdummy6to12

  subroutine dgdummy6to12(this, q12_t, q6_t, q12_0, q6_0, dg)
    
    implicit none

    class(constraint_6to12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q12_t, q12_0
    real(kind = 8), dimension( 6), intent(in) :: q6_t, q6_0

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)
    
  end subroutine dgdummy6to12

  subroutine kdummy6to12(this, lambda, k)

    implicit none

    class(constraint_6to12), intent(inout) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)

    real(kind = 8), allocatable, intent(inout) :: k(:,:)

  end subroutine kdummy6to12
  
  type(constraint_6to12_rigidtransition) function create_constraint_6to12_rigidtransition()
    
    create_constraint_6to12_rigidtransition%rankcount = 3
    create_constraint_6to12_rigidtransition%stiffness = .TRUE.
    return
  end function
  
  type(constraint_6to12_softtransition) function create_constraint_6to12_softtransition()
    
    create_constraint_6to12_softtransition%rankcount = 2
    create_constraint_6to12_softtransition%stiffness = .TRUE.
    return
  end function
  
  subroutine grigidtransition6to12(this, q12_t, q6_t, q12_0, q6_0, g)

    implicit none

    class(constraint_6to12_rigidtransition), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q12_t, q12_0
    real(kind = 8), dimension(6), intent(in) :: q6_t, q6_0
    
    real(kind = 8), allocatable, intent(inout) :: g(:)

    real(kind = 8), dimension(3) :: deltaphi_t, deltaphi_0

    deltaphi_0 = q6_0(1:3)-q12_0(1:3)
    deltaphi_t = q6_t(1:3)-q12_t(1:3)
       
    g(1) = dot_product(deltaphi_t, q12_t( 4: 6))-dot_product(deltaphi_0, q12_0( 4: 6))
    g(2) = dot_product(deltaphi_t, q12_t( 7: 9))-dot_product(deltaphi_0, q12_0( 7: 9))
    g(3) = dot_product(deltaphi_t, q12_t(10:12))-dot_product(deltaphi_0, q12_0(10:12))
    
  end subroutine grigidtransition6to12

  subroutine dgrigidtransition6to12(this, q12_t, q6_t, q12_0, q6_0, dg)

    implicit none

    class(constraint_6to12_rigidtransition), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q12_t, q12_0
    real(kind = 8), dimension( 6), intent(in) :: q6_t, q6_0

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)

    real(kind = 8), dimension(3) :: deltaphi_t

    deltaphi_t = q6_t(1:3)-q12_t(1:3)

    dg(:, :) = 0.0d0
    
    dg(1,  1: 3) =-q12_t( 4: 6)
    dg(1,  4: 6) = deltaphi_t
    dg(1, 13:15) = q12_t( 4: 6)

    dg(2,  1: 3) =-q12_t( 7: 9)
    dg(2,  7: 9) = deltaphi_t
    dg(2, 13:15) = q12_t( 7: 9)

    dg(3,  1: 3) =-q12_t(10:12)
    dg(3, 10:12) = deltaphi_t
    dg(3, 13:15) = q12_t(10:12)
    
    return
    
  end subroutine dgrigidtransition6to12

  subroutine krigidtransition6to12(this, lambda, k)

    implicit none

    class(constraint_6to12_rigidtransition), intent(inout) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)

    real(kind = 8), allocatable, intent(inout) :: k(:,:)
   
    k(:, :) = 0.0d0

    k( 4: 6,  1: 3) =-lambda(1)*i_33
    k( 1: 3,  4: 6) =-lambda(1)*i_33
    k( 4: 6, 13:15) = lambda(1)*i_33
    k(13:15,  4: 6) = lambda(1)*i_33

    k( 7: 9,  1: 3) =-lambda(2)*i_33
    k( 1: 3,  7: 9) =-lambda(2)*i_33
    k( 7: 9, 13:15) = lambda(2)*i_33
    k(13:15,  7: 9) = lambda(2)*i_33
    
    k(10:12,  1: 3) =-lambda(3)*i_33
    k( 1: 3, 10:12) =-lambda(3)*i_33
    k(10:12, 13:15) = lambda(3)*i_33
    k(13:15, 10:12) = lambda(3)*i_33
            
    return
    
  end subroutine krigidtransition6to12
  
  subroutine gsofttransition6to12(this, q12_t, q6_t, q12_0, q6_0, g)

    implicit none

    class(constraint_6to12_softtransition), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q12_t, q12_0
    real(kind = 8), dimension( 6), intent(in) :: q6_t, q6_0
    
    real(kind = 8), allocatable, intent(inout) :: g(:)

    real(kind = 8), dimension(3) :: deltaphi_t, tv_t, deltaphi_0
    
    deltaphi_0 = q6_0(1:3)-q12_0(1:3)

    this%dir(1) = dot_product(deltaphi_0, q12_0( 4: 6))
    this%dir(2) = dot_product(deltaphi_0, q12_0( 7: 9))
    this%dir(3) = dot_product(deltaphi_0, q12_0(10:12))    
    
    deltaphi_t = q6_t(1:3)-q12_t(1:3)
    
    tv_t = this%dir(2)*q12_t(4:6)-this%dir(1)*q12_t(7:9)
  
    g(1) = dot_product(deltaphi_t, tv_t)
    g(2) = dot_product(deltaphi_t, q12_t(10:12))-dot_product(deltaphi_0, q12_0(10:12))
      
    return
    
  end subroutine gsofttransition6to12

  subroutine dgsofttransition6to12(this, q12_t, q6_t, q12_0, q6_0, dg)

    implicit none

    class(constraint_6to12_softtransition), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q12_t, q12_0
    real(kind = 8), dimension( 6), intent(in) :: q6_t, q6_0

    real(kind = 8), allocatable, intent(inout) :: dg(:,:)

    real(kind = 8), dimension(3) :: deltaphi_t, tv_t, deltaphi_0

    deltaphi_0 = q6_0(1:3)-q12_0(1:3)
    
    this%dir(1) = dot_product(deltaphi_0, q12_0( 4: 6))
    this%dir(2) = dot_product(deltaphi_0, q12_0( 7: 9))    
    this%dir(3) = dot_product(deltaphi_0, q12_0(10:12))
    
    tv_t = this%dir(2)*q12_t(4:6)-this%dir(1)*q12_t(7:9)

    deltaphi_t = q6_t(1:3)-q12_t(1:3)

    dg(:, :) = 0.0d0
    
    dg(1,  1: 3) =-tv_t
    dg(1,  4: 6) = this%dir(2)*deltaphi_t
    dg(1,  7: 9) =-this%dir(1)*deltaphi_T
    dg(1, 13:15) = tv_t

    dg(2,  1: 3) =-q12_t(10:12)
    dg(2, 10:12) = deltaphi_t
    dg(2, 13:15) = q12_t(10:12)
        
    return
    
  end subroutine dgsofttransition6to12

  subroutine ksofttransition6to12(this, lambda, k)

    implicit none

    class(constraint_6to12_softtransition), intent(inout) :: this
    real(kind = 8), allocatable, intent(in) :: lambda(:)

    real(kind = 8), allocatable, intent(inout) :: k(:,:)
   
    k(:, :) = 0.0d0

    k( 1: 3,  4: 6) =-lambda(1)*this%dir(2)*i_33
    k( 4: 6,  1: 3) =-lambda(1)*this%dir(2)*i_33

    k( 1: 3,  7: 9) = lambda(1)*this%dir(1)*i_33
    k( 4: 6,  7: 9) = lambda(1)*this%dir(1)*i_33

    k( 4: 6, 13:15) = lambda(1)*this%dir(2)*i_33
    k(13:15,  4: 6) = lambda(1)*this%dir(2)*i_33

    k( 7: 9, 13:15) =-lambda(1)*this%dir(1)*i_33
    k(13:15,  7: 9) =-lambda(1)*this%dir(1)*i_33

    k( 1: 3, 10:12) =-lambda(2)*i_33
    k(10:12,  1: 3) =-lambda(2)*i_33

    k(10:12, 13:15) = lambda(2)*i_33
    k(13:15, 10:12) = lambda(2)*i_33
    
    return
    
  end subroutine ksofttransition6to12
      
end module class_constraint_6_to_12
