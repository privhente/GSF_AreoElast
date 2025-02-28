module my_materials
  
  use my_constants_structure, only: i_33
  implicit none

contains
  
  function isotropicmaterial(bige, nu) result(cscriptcart)
    
    implicit none
    
    real(kind = 8), intent(in) :: bige, nu
    
    real(kind = 8) :: cscriptcart(6, 6)
    
    cscriptcart(:, :) = 0.0d0
    
    cscriptcart(1, 1) = 1.0d0-nu
    cscriptcart(1, 2) = nu
    cscriptcart(1, 3) = nu
    
    cscriptcart(2, 1) = nu
    cscriptcart(2, 2) = 1.0d0-nu
    cscriptcart(2, 3) = nu
    
    cscriptcart(3, 1) = nu
    cscriptcart(3, 2) = nu
    cscriptcart(3, 3) = 1.0d0-nu
    
    cscriptcart(4:6, 4:6) = i_33*(1.0d0-2.0d0*nu)*0.5d0
    
    cscriptcart = bige/(1.0d0+nu)/(1.0d0-2.0d0*nu)*cscriptcart
    
    return
    
  end function isotropicmaterial

  function orthotropicmaterial(cmatv) result(cscriptcart)
    
    implicit none
    
    real(kind = 8), intent(in) :: cmatv(:)
    
    real(kind = 8) :: cscriptcart(6, 6)

    real(kind = 8) :: young1, young2, young3, shear23, shear13, shear12, poisson23, poisson13, poisson12

    real(kind = 8) :: work(6)
    integer :: ipiv(6), info

    cscriptcart(:, :) = 0.0d0

    young1    = cmatv(1)
    young2    = cmatv(2)
    young3    = cmatv(3)
    shear23   = cmatv(4)
    shear13   = cmatv(5)
    shear12   = cmatv(6)
    poisson23 = cmatv(7)
    poisson13 = cmatv(8)
    poisson12 = cmatv(9)
    
    cscriptcart(1, 1) = 1.0d0/young1
    cscriptcart(1, 2) = -poisson12/young1
    cscriptcart(1, 3) = -poisson13/young1
    
    cscriptcart(2, 1) = -poisson12/young1
    cscriptcart(2, 2) = 1.0d0/young2
    cscriptcart(2, 3) = -poisson23/young2
    
    cscriptcart(3, 1) = -poisson13/young1
    cscriptcart(3, 2) = -poisson23/young2
    cscriptcart(3, 3) = 1.0d0/young3

    cscriptcart(4, 4) = 1.0d0/shear23
    cscriptcart(5, 5) = 1.0d0/shear13
    cscriptcart(6, 6) = 1.0d0/shear12

!!$1000 format(1000f20.10)
!!$
!!$    print*, 'mat1'
!!$
!!$    do i = 1, 6
!!$
!!$       write(*, 1000) cscriptcart(i, :)
!!$       
!!$    end do

    call dgetrf(6, 6, cscriptcart, 6, ipiv, info)

    call dgetri(6, cscriptcart, 6, ipiv, work, 6, info)

!!$    print*, 'mat1'
!!$
!!$    do i = 1, 6
!!$
!!$       write(*, 1000) cscriptcart(i, :)
!!$       
!!$    end do
!!$
!!$pause
        
    return
    
  end function orthotropicmaterial
  
end module my_materials
