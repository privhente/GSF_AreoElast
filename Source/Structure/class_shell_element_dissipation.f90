module class_shell_element_dissipation

  implicit none

  type shell_element_dissipation

     real(kind = 8) :: e
     real(kind = 8) :: v
     real(kind = 8) :: energy
     real(kind = 8) :: sdiss(6)
     real(kind = 8) :: vdiss(3)
     real(kind = 8) :: cscript_ee(6, 6)
     real(kind = 8) :: cscript_vv(3, 3)
     real(kind = 8) :: cscript_ev(6, 3)
     real(kind = 8) :: cscript_ve(3, 6)

   contains

     procedure :: evaluate
     
  end type shell_element_dissipation
  
contains

  subroutine evaluate(this, e_1, e_2, v_1, v_2, cscript, rho, alpha_s, alpha_v, time, length_opt, tolerance_opt)

    use my_math_structure, only: outer, eye
  
    implicit none

    class(shell_element_dissipation) :: this
    
    real(kind = 8), dimension(6), intent(in) :: e_1, e_2
    real(kind = 8), dimension(3), intent(in) :: v_1, v_2
    real(kind = 8), dimension(6, 6), intent(in) :: cscript
    real(kind = 8), intent(in) :: rho, alpha_s, alpha_v, time
    real(kind = 8), optional, intent(in) :: length_opt
    real(kind = 8), optional, intent(in) :: tolerance_opt
    real(kind = 8) :: deltae(6), cscript_times_deltae(6), deltaes_cscript
    real(kind = 8) :: vs_1, vs_2, v_mean(3)
    real(kind = 8) :: length, tolerance

    length = 1.0d0
    tolerance = epsilon(1.0d0)

    if (present(length_opt   )) length    = length_opt
    if (present(tolerance_opt)) tolerance = tolerance_opt
    
    !!
    deltae = e_2-e_1
    cscript_times_deltae = matmul(cscript, deltae)
    deltaes_cscript = dsqrt(dot_product(deltae, cscript_times_deltae))
    
    !! 
    v_mean        = 0.5d0*(v_2+v_1)
    vs_1          = norm2(v_1)
    vs_2          = norm2(v_2)

    this%sdiss         = alpha_s/8.0d0*cscript_times_deltae
    this%cscript_ee    = alpha_s/8.0d0*cscript
        
    if (vs_2 > tolerance) then    
       this%vdiss(:)         = alpha_v/8.0d0* (vs_2-vs_1)/(vs_2+vs_1)*v_mean
       this%cscript_vv(:, :) = alpha_v/8.0d0*((vs_2-vs_1)/(vs_2+vs_1)*eye(3)*0.5d0+2.0d0*(vs_1/vs_2)/(vs_2+vs_1)**2*outer(v_mean, v_2))       
    else
       if (vs_1 > tolerance) then
          this%vdiss(:)         =-alpha_v/8.0d0*v_mean
          this%cscript_vv(:, :) =-alpha_v/8.0d0*eye(3)*0.5d0
       else
          this%vdiss(:)         = 0.0d0
          this%cscript_vv(:, :) = 0.0d0
       end if
    end if
    this%cscript_ev(:, :) = 0.0d0
    this%cscript_ve(:, :) = 0.0d0
    return
    
  end subroutine evaluate
  
end module class_shell_element_dissipation
