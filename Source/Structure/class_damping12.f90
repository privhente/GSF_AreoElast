module class_damping12
  use class_condensed12
  use my_constants_structure, only: i_33, e1, e2, e3
  use my_math_structure, only: outer, cross, skew
  
  implicit none

! ------------------------------------------------------------------------------  
 type, extends(condensed_12) :: damping12
    contains
      procedure :: damping12_f  => f_damping12_dummy
      procedure :: damping12_k  => k_damping12_dummy
  end type damping12
! ------------------------------------------------------------------------------  
  type damping12_c   !! container type                                               
    class(damping12), pointer :: c
  end type
! ------------------------------------------------------------------------------  
  type, extends(damping12) :: damping12_angularvelocity_local
  contains
    procedure :: damping12_f  => f_damping12_angularvelocity_local
    procedure :: damping12_k    => k_damping12_angularvelocity_local
  end type damping12_angularvelocity_local
  
  interface damping12_angularvelocity_local
    module procedure create_damping12_angularvelocity_local
  end interface
  
! ------------------------------------------------------------------------------  
  type, extends(damping12) :: damping12_angularvelocity_global
  contains
    procedure :: damping12_f  => f_damping12_angularvelocity_global
    procedure :: damping12_k    => k_damping12_angularvelocity_global
  end type damping12_angularvelocity_global
  
  interface damping12_angularvelocity_global
    module procedure create_damping12_angularvelocity_global
  end interface 
  
! ------------------------------------------------------------------------------  
  type, extends(damping12) :: damping12_translatory_local
  contains
    procedure :: damping12_f  => f_damping12_translatory_local
    procedure :: damping12_k    => k_damping12_translatory_local
  end type damping12_translatory_local

  interface damping12_translatory_local
    module procedure create_damping12_translatory_local
  end interface   
! ------------------------------------------------------------------------------  

  contains
  
! Include subroutines  
! ============================================================================
! ==== Dummy ==================================================================
! ============================================================================
  subroutine f_damping12_dummy(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    this%fqint(:) = 0.0d0
  end subroutine f_damping12_dummy
! ----------------------------------------------------------------------------
  subroutine k_damping12_dummy(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    this%kqq(:,:)     = 0.0d0
    this%kqv(:,:)     = 0.0d0
    this%penergy = 0.0d0
  end subroutine k_damping12_dummy
  
! ============================================================================
! ============================================================================  
! Functions for interfaces when using defined pointer
  type(damping12_angularvelocity_local) function create_damping12_angularvelocity_local()
    create_damping12_angularvelocity_local%node = 0
    return
  end function

  type(damping12_angularvelocity_global) function create_damping12_angularvelocity_global()
    create_damping12_angularvelocity_global%node = 0
    return
  end function

  type(damping12_translatory_local) function create_damping12_translatory_local()
    create_damping12_translatory_local%node = 0
    return
  end function
  
! ============================================================================
! ==== damping12_angularvelocity_local =======================================
! ============================================================================
  subroutine f_damping12_angularvelocity_local(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12_angularvelocity_local), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: w1(3), w2(3), w3(3), omega(3), Momega(3)
    real(kind = 8) :: d1(3), d2(3), d3(3), test
    real(kind = 8) :: ndir(3)

    d1    = 0.5d0*(qtn1(4:6)   + qtn(4:6))
    d2    = 0.5d0*(qtn1(7:9)   + qtn(7:9))
    d3    = 0.5d0*(qtn1(10:12) + qtn(10:12))
    
    w1    = 0.5d0*(vtn1(4:6)   + vtn(4:6)) 
    w2    = 0.5d0*(vtn1(7:9)   + vtn(7:9)) 
    w3    = 0.5d0*(vtn1(10:12) + vtn(10:12))
    
    ndir  = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3
    omega = 0.5d0* ( cross(d1,w1) + cross(d2,w2) + cross(d3,w3) )
    
    call this%linearstiffness(abs(dot_product(ndir,omega)))
    Momega = this%stiffness*matmul(outer(ndir,ndir),omega)
    
    this%fqint        = 0.0d0
    this%fqint(4:6)   = 0.5d0*cross(Momega,d1)
    this%fqint(7:9)   = 0.5d0*cross(Momega,d2)
    this%fqint(10:12) = 0.5d0*cross(Momega,d3)
    return
  end subroutine f_damping12_angularvelocity_local
! ----------------------------------------------------------------------------
  subroutine k_damping12_angularvelocity_local(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12_angularvelocity_local), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: w1(3), w2(3), w3(3), omega(3)
    real(kind = 8) :: d1(3), d2(3), d3(3)
    real(kind = 8) :: ndir(3), Comega(3,3), Momega(3)
    
    d1    = 0.5d0*(qtn1(4:6)   + qtn(4:6))
    d2    = 0.5d0*(qtn1(7:9)   + qtn(7:9))
    d3    = 0.5d0*(qtn1(10:12) + qtn(10:12))

    w1    = 0.5d0*(vtn1(4:6)   + vtn(4:6))
    w2    = 0.5d0*(vtn1(7:9)   + vtn(7:9))
    w3    = 0.5d0*(vtn1(10:12) + vtn(10:12))
    
    omega  = 0.5d0* ( cross(d1,w1) + cross(d2,w2) + cross(d3,w3) )
    ndir   = this%dir(1)*d1 + this%dir(2)*d2 + this%dir(3)*d3
    Comega = this%stiffness*outer(ndir,ndir)
    Momega = matmul(Comega,omega)

    this%kqq              = 0.0d0
    this%kqq(4:6,4:6)     = 0.5d0*0.5d0*( skew(Momega)  + 0.5d0*matmul(skew(d1),matmul(Comega,skew(w1))) - 0.5d0*this%stiffness*matmul(skew(d1),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(1))   )
    this%kqq(4:6,7:9)     = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d1),matmul(Comega,skew(w2))) - 0.5d0*this%stiffness*matmul(skew(d1),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(2))   )
    this%kqq(4:6,10:12)   = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d1),matmul(Comega,skew(w3))) - 0.5d0*this%stiffness*matmul(skew(d1),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(3))   )
    
    this%kqq(7:9,4:6)     = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d2),matmul(Comega,skew(w1))) - 0.5d0*this%stiffness*matmul(skew(d2),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(1))   )
    this%kqq(7:9,7:9)     = 0.5d0*0.5d0*( skew(Momega)  + 0.5d0*matmul(skew(d2),matmul(Comega,skew(w2))) - 0.5d0*this%stiffness*matmul(skew(d2),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(2))   )
    this%kqq(7:9,10:12)   = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d2),matmul(Comega,skew(w3))) - 0.5d0*this%stiffness*matmul(skew(d2),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(3))   )

    this%kqq(10:12,4:6)   = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d3),matmul(Comega,skew(w1))) - 0.5d0*this%stiffness*matmul(skew(d3),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(1))   )
    this%kqq(10:12,7:9)   = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d3),matmul(Comega,skew(w2))) - 0.5d0*this%stiffness*matmul(skew(d3),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(2))   )
    this%kqq(10:12,10:12) = 0.5d0*0.5d0*( skew(Momega)  + 0.5d0*matmul(skew(d3),matmul(Comega,skew(w3))) - 0.5d0*this%stiffness*matmul(skew(d3),(i_33*dot_product(omega,ndir) + outer(ndir,omega) )*this%dir(3))   )
    
    this%kqv              = 0.0d0
    this%kqv(4:6,4:6)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d1),matmul(Comega,skew(d1))) )
    this%kqv(4:6,7:9)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d1),matmul(Comega,skew(d2))) )
    this%kqv(4:6,10:12)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d1),matmul(Comega,skew(d3))) )
    
    this%kqv(7:9,4:6)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d2),matmul(Comega,skew(d1))) )
    this%kqv(7:9,7:9)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d2),matmul(Comega,skew(d2))) )
    this%kqv(7:9,10:12)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d2),matmul(Comega,skew(d3))) )
    
    this%kqv(10:12,4:6)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d3),matmul(Comega,skew(d1))) )
    this%kqv(10:12,7:9)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d3),matmul(Comega,skew(d2))) )
    this%kqv(10:12,10:12) = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d3),matmul(Comega,skew(d3))) )
    
    this%penergy = 0.5d0*dot_product(qtn1-q0,matmul(this%kqq, qtn1-q0))
    return
  end subroutine k_damping12_angularvelocity_local

! ============================================================================
! ==== damping12_angularvelocity_global =======================================
! ============================================================================
  subroutine f_damping12_angularvelocity_global(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12_angularvelocity_global), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: w1(3), w2(3), w3(3), omega(3), Momega(3)
    real(kind = 8) :: d1(3), d2(3), d3(3), test
    real(kind = 8) :: ndir(3)
    
    d1    = 0.5d0*(qtn1(4:6)   + qtn(4:6))
    d2    = 0.5d0*(qtn1(7:9)   + qtn(7:9))
    d3    = 0.5d0*(qtn1(10:12) + qtn(10:12))
    
    w1    = 0.5d0*(vtn1(4:6)   + vtn(4:6)) 
    w2    = 0.5d0*(vtn1(7:9)   + vtn(7:9)) 
    w3    = 0.5d0*(vtn1(10:12) + vtn(10:12))
    
    ndir  = this%dir
    omega = 0.5d0* ( cross(d1,w1) + cross(d2,w2) + cross(d3,w3) )
    
    call this%linearstiffness(abs(dot_product(ndir,omega)))
    Momega = this%stiffness*matmul(outer(ndir,ndir),omega)
    
    this%fqint        = 0.0d0
    this%fqint(4:6)   = 0.5d0*cross(Momega,d1)
    this%fqint(7:9)   = 0.5d0*cross(Momega,d2)
    this%fqint(10:12) = 0.5d0*cross(Momega,d3)
    return
  end subroutine f_damping12_angularvelocity_global
! ----------------------------------------------------------------------------
  subroutine k_damping12_angularvelocity_global(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12_angularvelocity_global), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2
    real(kind = 8) :: w1(3), w2(3), w3(3), omega(3)
    real(kind = 8) :: d1(3), d2(3), d3(3)
    real(kind = 8) :: ndir(3), Comega(3,3), Momega(3)
    
    d1    = 0.5d0*(qtn1(4:6)   + qtn(4:6))
    d2    = 0.5d0*(qtn1(7:9)   + qtn(7:9))
    d3    = 0.5d0*(qtn1(10:12) + qtn(10:12))

    w1    = 0.5d0*(vtn1(4:6)   + vtn(4:6))
    w2    = 0.5d0*(vtn1(7:9)   + vtn(7:9))
    w3    = 0.5d0*(vtn1(10:12) + vtn(10:12))
    
    omega  = 0.5d0* ( cross(d1,w1) + cross(d2,w2) + cross(d3,w3) )
    ndir   = this%dir
    Comega = this%stiffness*outer(ndir,ndir)
    Momega = matmul(Comega,omega)

    this%kqq              = 0.0d0
    this%kqq(4:6,4:6)     = 0.5d0*0.5d0*( skew(Momega)  + 0.5d0*matmul(skew(d1),matmul(Comega,skew(w1))) )
    this%kqq(4:6,7:9)     = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d1),matmul(Comega,skew(w2))) )
    this%kqq(4:6,10:12)   = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d1),matmul(Comega,skew(w3))) )
    
    this%kqq(7:9,4:6)     = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d2),matmul(Comega,skew(w1))) )
    this%kqq(7:9,7:9)     = 0.5d0*0.5d0*( skew(Momega)  + 0.5d0*matmul(skew(d2),matmul(Comega,skew(w2))) )
    this%kqq(7:9,10:12)   = 0.5d0*0.5d0*(               - 0.5d0*matmul(skew(d2),matmul(Comega,skew(w3))) )

    this%kqq(10:12,4:6)   = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d3),matmul(Comega,skew(w1))) )
    this%kqq(10:12,7:9)   = 0.5d0*0.5d0*(               + 0.5d0*matmul(skew(d3),matmul(Comega,skew(w2))) )
    this%kqq(10:12,10:12) = 0.5d0*0.5d0*( skew(Momega)  + 0.5d0*matmul(skew(d3),matmul(Comega,skew(w3))) )
    
    this%kqv              = 0.0d0
    this%kqv(4:6,4:6)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d1),matmul(Comega,skew(d1))) )
    this%kqv(4:6,7:9)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d1),matmul(Comega,skew(d2))) )
    this%kqv(4:6,10:12)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d1),matmul(Comega,skew(d3))) )
    
    this%kqv(7:9,4:6)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d2),matmul(Comega,skew(d1))) )
    this%kqv(7:9,7:9)     = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d2),matmul(Comega,skew(d2))) )
    this%kqv(7:9,10:12)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d2),matmul(Comega,skew(d3))) )
    
    this%kqv(10:12,4:6)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d3),matmul(Comega,skew(d1))) )
    this%kqv(10:12,7:9)   = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d3),matmul(Comega,skew(d2))) )
    this%kqv(10:12,10:12) = 0.5d0*0.5d0*( - 0.5d0*matmul(skew(d3),matmul(Comega,skew(d3))) )
    
    this%penergy = 0.5d0*dot_product(qtn1-q0,matmul(this%kqq, qtn1-q0))
    return
  end subroutine k_damping12_angularvelocity_global
  
! ============================================================================
! ==== damping12_translatory_local ============================================
! ============================================================================
  subroutine f_damping12_translatory_local(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12_translatory_local), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2

    return
  end subroutine f_damping12_translatory_local
! ----------------------------------------------------------------------------
  subroutine k_damping12_translatory_local(this, vtn1, vtn, qtn1, qtn, q0, qopt1, qopt2)
    implicit none
    class(damping12_translatory_local), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: vtn1, vtn, qtn1, qtn, q0
    real(kind = 8), dimension(12), intent(in), optional :: qopt1, qopt2

    return
  end subroutine k_damping12_translatory_local
! ============================================================================  
end module class_damping12