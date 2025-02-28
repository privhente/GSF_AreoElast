module class_node_12
  
  use class_node_structure

  implicit none
  
!> @brief Node with three directors
!
!> position vector \f$\phi\f$ and directors \f$d_n\f$, using three spatial degrees of freedom
!> \f[
!>  q= \left( \phi \, d_1 \, d_2 \, d_3 \right)^T
!> \f]
  type, extends(node_structure) :: node_12
     
     integer :: coordinates(12), velocity(12)
     real(kind = 8) :: q_0(12)
     real(kind = 8) :: fqext(12)
     real(kind = 8) :: kextq(12, 12)
     real(kind = 8) :: ml_t(6)
   contains
     
     procedure :: externalterms => nodalexternalterms12
     procedure :: materialterms => nodalmaterialterms12
     
  end type node_12
  
contains
  
!> @brief calculate external nodal forces and moments
!
!> Calculates nodal forces and moments from spatial and material external forces and moments.
!> Also calculates the nodal stiffness matrix Kqq resulting from moment loading
!
!> @param[in] q_t Geometry configuration at current time (position and three directors
!> @param[in] ml_t material load vector at current time (components 1:3 are forces, 4:6 are moments)
!> @param[in] sl_t spatial load vector at current time (components 1:3 are forces, 4:6 are moments)
  
  subroutine nodalexternalterms12(this, q_t, ml_t, sl_t)
    
    implicit none

    class(node_12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q_t    
    real(kind = 8), dimension(6) , intent(in) :: ml_t
    real(kind = 8), dimension(6) , intent(in) :: sl_t
    
    real(kind = 8), dimension(3) :: d1_t, d2_t, d3_t
    real(kind = 8), dimension(3, 3) :: rotation
    real(kind = 8), dimension(3) :: mscript
    
    d1_t = q_t( 4: 6)
    d2_t = q_t( 7: 9)
    d3_t = q_t(10:12)       
   
    rotation = outer(d1_t, e1)+outer(d2_t, e2)+outer(d3_t, e3)
   
    mscript = sl_t(4:6)+matmul(rotation, ml_t(4:6))
    this%fqext( 1: 3) = sl_t(1:3) + matmul(rotation, ml_t(1:3))
    this%fqext( 4: 6) = 0.5d0*cross(mscript, d1_t)
    this%fqext( 7: 9) = 0.5d0*cross(mscript, d2_t)
    this%fqext(10:12) = 0.5d0*cross(mscript, d3_t)
              
    this%kextq(:, :) = 0.0d0
    this%kextq( 1: 3,  4: 6) = ml_t(1)*i_33
    this%kextq( 1: 3,  7: 9) = ml_t(2)*i_33
    this%kextq( 1: 3, 10:12) = ml_t(3)*i_33

    this%kextq( 4: 6,  4: 6) =-0.5d0*(-skew(mscript) + ml_t(4)*skew(d1_t))
    this%kextq( 4: 6,  7: 9) =-0.5d0*                ml_t(5)*skew(d1_t)
    this%kextq( 4: 6, 10:12) =-0.5d0*                ml_t(6)*skew(d1_t)

    this%kextq( 7: 9,  4: 6) =-0.5d0*                ml_t(4)*skew(d2_t)
    this%kextq( 7: 9,  7: 9) =-0.5d0*(-skew(mscript) + ml_t(5)*skew(d2_t))
    this%kextq( 7: 9, 10:12) =-0.5d0*                ml_t(6)*skew(d2_t)

    this%kextq(10:12,  4: 6) =-0.5d0*                ml_t(4)*skew(d3_t)
    this%kextq(10:12,  7: 9) =-0.5d0*                ml_t(5)*skew(d3_t)
    this%kextq(10:12, 10:12) =-0.5d0*(-skew(mscript) + ml_t(6)*skew(d3_t))
   
    return

  end subroutine nodalexternalterms12

!> @brief convert the external nodal forces and moments in director-based components to material forces and moments in rotation-based terms
!
!> @param[in] q_t Geometry configuration at current time (position and three directors
!> @param[in] fqext external load vector in director-based terms
!> @param[in] ml_t material load vector at current time (components 1:3 are forces, 4:6 are moments)
  
  subroutine nodalmaterialterms12(this, q_t, fqext)
    
    implicit none

    class(node_12), intent(inout) :: this
    real(kind = 8), dimension(12), intent(in) :: q_t, fqext  
    
    real(kind = 8), dimension(3)    :: d1_t, d2_t, d3_t
    real(kind = 8), dimension(3)    :: f_x, f_d1, f_d2, f_d3
    real(kind = 8), dimension(3, 3) :: rotation
    
    d1_t = q_t( 4: 6)
    d2_t = q_t( 7: 9)
    d3_t = q_t(10:12)            
    rotation = outer(d1_t, e1)+outer(d2_t, e2)+outer(d3_t, e3)

    !< determining the material load terms
    f_x  = fqext(1:3)
    f_d1 = fqext(4:6)
    f_d2 = fqext(7:9)
    f_d3 = fqext(10:12)
      
    this%ml_t(1:3) = matmul(transpose(rotation),f_x)
    this%ml_t(4:6) = matmul(transpose(rotation),cross(f_d1,d1_t)+cross(f_d2,d2_t)+cross(f_d3,d3_t))
   
    return

  end subroutine nodalmaterialterms12
  
end module class_node_12
