module class_spring12
  use class_condensed12
  use my_constants_structure, only: i_33, e1, e2, e3
  use my_math_structure, only: outer
  
  implicit none

! ------------------------------------------------------------------------------  
 type, extends(condensed_12) :: spring12
    contains
      procedure :: spring12_f  => f_spring12_dummy
      procedure :: spring12_k    => k_spring12_dummy
  end type spring12
! ------------------------------------------------------------------------------  
  type spring12_c   !! container type                                               
    class(spring12), pointer :: c
  end type
! ------------------------------------------------------------------------------  
  type, extends(spring12) :: spring12_translatory_inplane
  contains
    procedure :: spring12_f  => f_spring12_translatory_inplane
    procedure :: spring12_k    => k_spring12_translatory_inplane
  end type spring12_translatory_inplane
  
  interface spring12_translatory_inplane
    module procedure create_spring12_translatory_inplane
  end interface  
! ------------------------------------------------------------------------------  
  type, extends(spring12) :: spring12_translatory_global
  contains
    procedure :: spring12_f  => f_spring12_translatory_global
    procedure :: spring12_k    => k_spring12_translatory_global
  end type spring12_translatory_global
  
  interface spring12_translatory_global
    module procedure create_spring12_translatory_global
  end interface    
! ------------------------------------------------------------------------------  
  type, extends(spring12) :: spring12_translatory_local
  contains
    procedure :: spring12_f  => f_spring12_translatory_local
    procedure :: spring12_k    => k_spring12_translatory_local
  end type spring12_translatory_local
  
  interface spring12_translatory_local
    module procedure create_spring12_translatory_local
  end interface   
! ------------------------------------------------------------------------------  
  type, extends(spring12) :: spring12_translatory_corotational
  contains
    procedure :: spring12_f  => f_spring12_translatory_corotational
    procedure :: spring12_k    => k_spring12_translatory_corotational
  end type spring12_translatory_corotational

  interface spring12_translatory_corotational
    module procedure create_spring12_translatory_corotational
  end interface     
  contains
  
! Include subroutines  
! ============================================================================
! ==== Dummy ==================================================================
! ============================================================================
  subroutine f_spring12_dummy(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    this%fqint(:) = 0.0d0
  end subroutine f_spring12_dummy
! ----------------------------------------------------------------------------
  subroutine k_spring12_dummy(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    this%kqq(:,:) = 0.0d0
    this%penergy  = 0.0d0
  end subroutine k_spring12_dummy
  
! ============================================================================
! ============================================================================  
! Functions for interfaces when using defined pointer
  type(spring12_translatory_inplane) function create_spring12_translatory_inplane()
    create_spring12_translatory_inplane%node = 0
    return
  end function

  type(spring12_translatory_global) function create_spring12_translatory_global()
    create_spring12_translatory_global%node = 0
    return
  end function

  type(spring12_translatory_local) function create_spring12_translatory_local()
    create_spring12_translatory_local%node = 0
    return
  end function
  
  type(spring12_translatory_corotational) function create_spring12_translatory_corotational()
    create_spring12_translatory_corotational%node = 0
    return
  end function
  
! ============================================================================
! ==== spring12_translatory_global ===========================================
! ============================================================================
  subroutine f_spring12_translatory_global(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_global), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: delta_x(3), ndir(3)
    this%fqint = 0.0d0
    delta_x(:) = 0.5d0*(qtn1(1:3) + qtn(1:3)) - q0(1:3)
    ndir = this%dir

    call this%linearstiffness(abs(dot_product(ndir,delta_x)))
    this%fqint(1:3) = this%stiffness*matmul(outer(ndir,ndir),delta_x)
    return
  end subroutine f_spring12_translatory_global
! ----------------------------------------------------------------------------
  subroutine k_spring12_translatory_global(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_global), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: ndir(3)
    this%kqq = 0.0d0
    ndir = this%dir
    
    this%kqq(1:3,1:3) = 0.5d0*this%stiffness*outer(ndir,ndir)
    this%penergy = 0.5d0*dot_product(qtn1(1:3)-q0(1:3),matmul(this%kqq(1:3,1:3), qtn1(1:3)-q0(1:3)))
    return
  end subroutine k_spring12_translatory_global
! ============================================================================
! ==== spring12_translatory_local ============================================
! ============================================================================
  subroutine f_spring12_translatory_local(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_local), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: delta_x(3), ndir(3), d1(3), d2(3), d3(3)
    this%fqint = 0.0d0
    delta_x(:) = 0.5d0*(qtn1(1:3)   + qtn(1:3)) - q0(1:3)
    d1         = 0.5d0*(qtn1(4:6)   + qtn(4:6))
    d2         = 0.5d0*(qtn1(7:9)   + qtn(7:9))
    d3         = 0.5d0*(qtn1(10:12) + qtn(10:12))
    ndir       = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3
    
    call this%linearstiffness(abs(dot_product(ndir,delta_x)))
    this%fqint(1:3) = this%stiffness* (matmul(outer(ndir,ndir),delta_x))
    return
  end subroutine f_spring12_translatory_local
! ----------------------------------------------------------------------------
  subroutine k_spring12_translatory_local(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_local), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: delta_x(3), ndir(3), d1(3), d2(3), d3(3)
    real(kind = 8) :: A_n(3,9)
    this%kqq = 0.0d0
    delta_x(:) = 0.5d0*(qtn1(1:3)   + qtn(1:3)) - q0(1:3)
    d1         = 0.5d0*(qtn1(4:6)   + qtn(4:6))
    d2         = 0.5d0*(qtn1(7:9)   + qtn(7:9))
    d3         = 0.5d0*(qtn1(10:12) + qtn(10:12))
    ndir       = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3

    A_n(:,:)     = 0.0d0
    A_n(1:3,1:3) = i_33*this%dir(1)
    A_n(1:3,4:6) = i_33*this%dir(2)
    A_n(1:3,7:9) = i_33*this%dir(3)
    
    this%kqq(:,:)      = 0.0d0
    this%kqq(1:3,1:3)  = 0.5d0*this%stiffness*outer(ndir,ndir)
    this%kqq(1:3,4:12) = 0.5d0*this%stiffness*(matmul(outer(ndir,delta_x),A_n) + dot_product(ndir,delta_x)*A_n)
    this%penergy  = 0.5d0*dot_product(qtn1-q0,matmul(this%kqq, qtn1-q0))    
    return
  end subroutine k_spring12_translatory_local
! ============================================================================  
! ==== spring12_translatory_corotational =====================================
! ============================================================================
  subroutine f_spring12_translatory_corotational(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_corotational), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: delta_x(3), ndir(3)
    this%fqint = 0.0d0
    delta_x(:) = 0.5d0*(qtn1(1:3)   + qtn(1:3)) - q0(1:3)
    if (norm2(delta_x) .eq. 0.0d0) then
      return
    end if

    ndir = delta_x/norm2(delta_x)
    
    call this%linearstiffness(abs(dot_product(ndir,delta_x)))
    this%fqint(1:3) = this%stiffness* (matmul(outer(ndir,ndir),delta_x))
    return
  end subroutine f_spring12_translatory_corotational
! ----------------------------------------------------------------------------
  subroutine k_spring12_translatory_corotational(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_corotational), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: delta_x(3), ndir(3), A_n(3,3)
    this%kqq = 0.0d0
    delta_x(:) = 0.5d0*(qtn1(1:3) + qtn(1:3)) - q0(1:3)

    if (norm2(delta_x) .eq. 0.0d0) then
      return
    end if

    ndir     = delta_x/norm2(delta_x)    
    A_n(:,:) = i_33/norm2(delta_x) - outer(delta_x,delta_x)/(dot_product(delta_x,delta_x)**(1.5d0))
    
    this%kqq(1:3,1:3) = 0.5d0*this%stiffness*(matmul(outer(ndir,delta_x),A_n) + outer(ndir,ndir) + dot_product(ndir,delta_x)*A_n)
    this%penergy = 0.5d0*dot_product(qtn1(1:3)-q0(1:3),matmul(this%kqq(1:3,1:3), qtn1(1:3)-q0(1:3)))
    return
  end subroutine k_spring12_translatory_corotational
! ============================================================================    
! ==== spring12_translatory_inplane_ =========================================
! ============================================================================
  subroutine f_spring12_translatory_inplane(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_inplane), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: ndir(3), nxdir(3), delta_x(3)
    this%fqint = 0.0d0
    delta_x(:) = 0.5d0*(qtn1(1:3) + qtn(1:3)) - q0(1:3)

    nxdir(1) = this%dir(1)*delta_x(1)
    nxdir(2) = this%dir(2)*delta_x(2)
    nxdir(3) = this%dir(3)*delta_x(3)
    
    if (norm2(nxdir) .eq. 0.0d0) then
      ndir = this%dir
    else
      ndir = nxdir/norm2(nxdir)
    end if
    
    call this%linearstiffness(abs(dot_product(ndir,delta_x)))
    this%fqint(1:3) = this%stiffness*matmul(outer(ndir,ndir),delta_x)
    return
  end subroutine f_spring12_translatory_inplane
! ----------------------------------------------------------------------------
  subroutine k_spring12_translatory_inplane(this, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(spring12_translatory_inplane), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: ndir(3), nxdir(3), delta_x(3), N_x(3,3), delta_N(3,3)
    this%kqq = 0.0d0
    delta_x(:) = 0.5d0*(qtn1(1:3) + qtn(1:3)) - q0(1:3)

    nxdir(1) = this%dir(1)*delta_x(1)
    nxdir(2) = this%dir(2)*delta_x(2)
    nxdir(3) = this%dir(3)*delta_x(3)
    delta_N  = 0.0d0
    N_x      = 0.0d0
    
    if (norm2(nxdir) .eq. 0.0d0) then
      ndir = this%dir
    else
      ndir = nxdir/norm2(nxdir)
      N_x(1,1)     = this%dir(1)
      N_x(2,2)     = this%dir(2)
      N_x(3,3)     = this%dir(3)
      delta_N(:,:) = matmul(i_33/norm2(nxdir) - outer(nxdir,nxdir)/(dot_product(nxdir,nxdir)**(1.5d0)),N_x)
    end if
    
    call this%linearstiffness(abs(dot_product(ndir,delta_x)))
    this%kqq(1:3,1:3) = 0.5d0*this%stiffness*(matmul(outer(ndir,delta_x),delta_N) + outer(ndir,ndir) + dot_product(ndir,delta_x)*delta_N)
    this%penergy = 0.5d0*dot_product(qtn1(1:3)-q0(1:3),matmul(this%kqq(1:3,1:3), qtn1(1:3)-q0(1:3)))
    return
  end subroutine k_spring12_translatory_inplane
! ============================================================================
end module class_spring12