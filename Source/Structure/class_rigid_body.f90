module class_rigid_body

  use my_constants_structure, only: i_33
  use my_math_structure, only: eye
  use class_node_12
  
  implicit none
  
  type :: rigid_body

     type(node_12) :: node
     integer :: property
     integer :: indicesq12(12)
     integer :: nnodes = 1
     real(kind = 8) :: cmass(4, 4)
     real(kind = 8), dimension(12) :: f, fqint, fqdyn, fv, fqext
     real(kind = 8), dimension(12, 12) :: m, mhat, kqq, kvv, kqv, kvq, keqq
     real(kind = 8) :: q_0(12)
     real(kind = 8) :: deltat
     integer :: simutype
     real(kind = 8) :: alpha
     real(kind = 8) :: penergy
     real(kind = 8) :: denergy
     logical :: flag_kgeo_on, flag_kmat_on

  contains

     procedure :: initialization => rigidbodyinitialization
     procedure :: actualization  => rigidbodyactualization     
     procedure, private :: rigidbodym
     procedure, private :: rigidbodyfk
     
  end type rigid_body

contains

!!$! This subroutine constructs the object and evaluates m.

  subroutine rigidbodyinitialization(this)

    implicit none
    class(rigid_body), intent(inout) :: this
    call this%rigidbodym()
    return
    
  end subroutine rigidbodyinitialization

! This subroutine evaluates f and k for q_t, mload and sload.
  
  subroutine rigidbodyactualization(this, q_1, ml_a, sl_a, q_2, v_1, v_2)

    implicit none
    class(rigid_body) , intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q_1(12)
    real(kind = 8), dimension(6), intent(in) :: ml_a
    real(kind = 8), dimension(6), intent(in) :: sl_a 
    real(kind = 8), dimension(12), optional, intent(in) :: q_2(12)
    real(kind = 8), dimension(12), optional, intent(in) :: v_1(12)
    real(kind = 8), dimension(12), optional, intent(in) :: v_2(12)
    
    if (present(q_2)) then
       if (present(v_1) .and. present(v_2)) then
          call this%rigidbodyfk(q_1, ml_a, sl_a, q_2, v_1, v_2)
       else
          call this%rigidbodyfk(q_1, ml_a, sl_a, q_2)
       end if
    else
       call this%rigidbodyfk(q_1, ml_a, sl_a)
    end if
    
    return
    
  end subroutine rigidbodyactualization
    
! This function computes the mass matrix.
  subroutine rigidbodym(this)

    implicit none
    
    class(rigid_body), intent(inout) :: this
    real(kind = 8), dimension(12, 12) :: m
    
    m( 1: 3,  1: 3) = this%cmass(1, 1)*i_33
    m( 1: 3,  4: 6) = this%cmass(1, 2)*i_33
    m( 1: 3,  7: 9) = this%cmass(1, 3)*i_33
    m( 1: 3, 10:12) = this%cmass(1, 4)*i_33
    
    m( 4: 6,  1: 3) = this%cmass(2, 1)*i_33
    m( 4: 6,  4: 6) = this%cmass(2, 2)*i_33
    m( 4: 6,  7: 9) = this%cmass(2, 3)*i_33
    m( 4: 6, 10:12) = this%cmass(2, 4)*i_33
    
    m( 7: 9,  1: 3) = this%cmass(3, 1)*i_33
    m( 7: 9,  4: 6) = this%cmass(3, 2)*i_33
    m( 7: 9,  7: 9) = this%cmass(3, 3)*i_33
    m( 7: 9, 10:12) = this%cmass(3, 4)*i_33
    
    m(10:12,  1: 3) = this%cmass(4, 1)*i_33
    m(10:12,  4: 6) = this%cmass(4, 2)*i_33
    m(10:12,  7: 9) = this%cmass(4, 3)*i_33
    m(10:12, 10:12) = this%cmass(4, 4)*i_33

    this%m = m
   
    return

  end subroutine rigidbodym

! This subroutine computes the generalized force vector due to external loads and its corresponding stiffness matrix. 
  subroutine rigidbodyfk(this, q_1, ml_a, sl_a, q_2_opt, v_1_opt, v_2_opt)

    implicit none
    class(rigid_body), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q_1
    real(kind = 8), dimension(6), intent(in) :: ml_a, sl_a
    real(kind = 8), dimension(12), optional, intent(in) :: q_2_opt
    real(kind = 8), dimension(12), optional, intent(in) :: v_1_opt
    real(kind = 8), dimension(12), optional, intent(in) :: v_2_opt
    real(kind = 8), dimension(12) :: q_2, v_1, v_2
    
    if (present(q_2_opt)) then
       q_2 = q_2_opt
    else
       q_2 = q_1
    end if
        
    if (present(v_1_opt) .and. present(v_2_opt)) then
       v_1 = v_1_opt
       v_2 = v_2_opt
    else
       v_1(:) = 0.0d0
       v_2(:) = 0.0d0
    end if
    
    this%penergy   = 0.0d0
    this%fqint(:)  = 0.0d0
    this%fqdyn(:)  = 0.0d0
    this%fv(:)     = 0.0d0    
    this%kqq(:, :) = 0.0d0
    this%kvv(:, :) = 0.0d0
    this%kqv(:, :) = 0.0d0
    this%kvq(:, :) = 0.0d0

    this%fqdyn   = matmul(this%m, (v_2-v_1)/this%deltat)
    this%fqint   = 0.0d0
    this%fv      = matmul(this%m, (q_2-q_1)/this%deltat-0.5d0*(v_1+v_2))
    
    this%kqv     = this%m/this%deltat
    this%kvv     =-this%m*0.5d0
    this%kvq     = this%m/this%deltat
    
! evaluation of external terms, i.e. material and spatial loads.
    if (present(q_2_opt)) then
       call this%node%externalterms(0.5d0*(q_1+q_2), ml_a, sl_a)
    else
       call this%node%externalterms(q_1, ml_a, sl_a)
    end if
    
    this%fqext = this%node%fqext
    this%kqq   = this%kqq + 0.5d0*this%node%kextq
    
    return
    
  end subroutine rigidbodyfk
    
end module class_rigid_body
