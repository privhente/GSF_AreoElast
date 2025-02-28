module class_beam_element_dissipation

  implicit none

  type beam_element_dissipation
    real(kind = 8) :: s_diss(6)
    real(kind = 8) :: C_diss(6,6)
    real(kind = 8) :: v_diss(24)
    real(kind = 8) :: kvv_diss(24,24)

    contains
      procedure :: evaluate
  end type beam_element_dissipation
  
contains

  subroutine evaluate(this, strain_1, strain_2, v_1, v_2, c_beam, m, alpha_s, alpha_v)
    use my_math_structure, only: outer, eye
    implicit none
    class(beam_element_dissipation), intent(inout):: this
    
    real(kind = 8), dimension(6), intent(in) :: strain_1, strain_2
    real(kind = 8), dimension(24), intent(in) :: v_1, v_2
    real(kind = 8), dimension(6, 6), intent(in) :: c_beam
    real(kind = 8), intent(in) :: alpha_s, alpha_v
    real(kind = 8), dimension(24,24), intent(in) :: m       ! mass matrix
    
    real(kind = 8) :: Dv, hv, Dv_v(24), hv_v(24)            ! energy function and denominator for vdiss
    real(kind = 8) :: normM_v_1, normM_v_2                  ! Norm of v*m*v
    real(kind = 8) :: normC_s                               ! E*C*E, Dissipation function for sdiss, denominator
    
    Dv = 0.0d0
    hv = 0.0d0
    Dv_v(:) = 0.0d0
    hv_v(:) = 0.0d0
    normM_v_1 = 0.0d0
    normM_v_2 = 0.0d0
    normC_s = 0.0d0
    
    ! Energy-Dissipation: aus Energy-dissipative momentum-conserving time-stepping algorithms for the dynamics of nonlinear cosserat rods (Armero+Romero 2003)
    ! Stress dissipation using 1st order scheme:
    if (alpha_s .ne. 0.0d0) then
      normC_s     = dsqrt(dot_product(strain_2-strain_1,matmul(c_beam,strain_2-strain_1)))
      this%s_diss = 0.5d0*alpha_s*matmul(c_beam,strain_2-strain_1)
      this%C_diss = 0.5d0*alpha_s*c_beam
    end if
    
    ! velocity dissipation using 1st order scheme: dissipation function = 1/2*alpha*( ||v_2||_m - ||v_1||_m )^2
    if (alpha_v .ne. 0.0d0) then
      normM_v_1 = dsqrt(dot_product(v_1,matmul(m,v_1)))
      normM_v_2 = dsqrt(dot_product(v_2,matmul(m,v_2)))
      Dv = (normM_v_2 - normM_v_1);
      hv = (normM_v_2 + normM_v_1);
      if ((hv .eq. 0.0d0) .or. (normM_v_2 .eq. 0.0d0)) then
        this%v_diss(:) = 0.0d0
        this%kvv_diss  = 0.0d0
        return
      end if
      this%v_diss(:) = 0.5d0*alpha_v*(v_2+v_1)*Dv/hv
      Dv_v = matmul(m,v_2)/normM_v_2
      hv_v = Dv_v
      this%kvv_diss  = 0.5d0*alpha_v*Dv/hv*eye(24) + 0.5d0*alpha_v*outer(v_2+v_1,(Dv_v*hv - Dv*hv_v)/(hv*hv))
    end if
    
    return
  end subroutine evaluate
  
end module class_beam_element_dissipation
