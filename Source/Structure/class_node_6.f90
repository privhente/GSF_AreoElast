module class_node_6
  
  use class_node_structure
  
  implicit none

  type, extends(node_structure) :: node_6
     
     integer :: coordinates(6), coordinates_q(6)
     real(kind = 8) :: q_0(6)
     real(kind = 8) :: fqext(6)
     real(kind = 8) :: kextq(6,6)
     
   contains

     procedure :: externalterms => nodalexternalterms6
     
  end type node_6
  
contains

!> @brief calculate external nodal forces and moments
!
!> Calculates nodal forces and moments from spatial and material external forces and moments.
!> Also calculates the nodal stiffness matrix Kqq resulting from moment loading
!
!> @param[in] q_t Geometry configuration at current time (position and three directors
!> @param[in] ml_t material load vector at current time (components 1:3 are forces, 4:6 are moments)
!> @param[in] sl_t spatial load vector at current time (components 1:3 are forces, 4:6 are moments)

  subroutine nodalexternalterms6(this, q_t, ml_t, sl_t)
    
    implicit none

    class(node_6), intent(inout) :: this
    real(kind = 8), dimension(6), intent(in) :: q_t
    real(kind = 8), dimension(6), intent(in) :: ml_t
    real(kind = 8), dimension(6), intent(in) :: sl_t
    
    real(kind = 8), dimension(3) :: dir_t

    dir_t = q_t( 4: 6)
    
    this%fqext(:) = 0.0d0

    this%kextq(:, :) = 0.0d0

    this%fqext(1:3) = sl_t(1:3)+ml_t(3)*dir_t

    this%kextq(1:3, 4:6) = i_33*ml_t(3)

    return

  end subroutine nodalexternalterms6

end module class_node_6
