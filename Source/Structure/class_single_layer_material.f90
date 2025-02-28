module class_single_layer_material

  use my_constants_structure, only: deg2rad, e1, e2, e3
  
  use my_materials, only: orthotropicmaterial
  
  use my_math_structure, only: curva2curvb3matrix, curva2curvb6matrix
  
  implicit none

  type :: single_layer_material ! single layer material for the extensible-director-based solid-degenerate shell.

     integer :: material ! material number.
     real(kind = 8) :: cscriptmat(6, 6) ! three-dimensional constitutive law in Voigth form expressed in the material coordinate system.
     real(kind = 8) :: cscript(6, 6) ! three-dimensional constitutive law in Voigth form expressed in the elemental coordinate system.
     real(kind = 8) :: thickness ! thickness of the layer.
     real(kind = 8) :: density ! mass densitiy per unit volume.
     real(kind = 8) :: alpha_v ! coefficient of vel. dissipation.
     real(kind = 8) :: alpha_s ! coefficient of stress dissipation.
     real(kind = 8) :: angledeg ! angle of orientation about the normal direction in degrees.
     real(kind = 8) :: angle ! angle of orientation about the normal direction in rad.
     real(kind = 8), dimension(3) :: e1mat, e2mat, e3mat ! basis of the material coordinate system.
     
   contains

     procedure :: initialization => initializationsinglelayermaterial
     
  end type single_layer_material

contains

  subroutine initializationsinglelayermaterial(this, cmatv, density, thickness,  angledeg_opt, alpha_opt_s, alpha_opt_v) ! this subroutine initializes the sinle layer material.    

    implicit none

    class(single_layer_material), intent(inout) :: this
    real(kind = 8), intent(in) :: cmatv(:) ! this vector contains E1, E2, E3, G23, G13, G12, nu23, nu13, nu12 (orthotropic material)
    real(kind = 8), intent(in) :: density ! this is the density per unit volume.
    real(kind = 8), intent(in) :: thickness ! this is the thickness of the layer.
    real(kind = 8), optional, intent(in) :: angledeg_opt ! angle of orientation about the normal direction in degrees. If not specified give just null().
    real(kind = 8), optional, intent(in) :: alpha_opt_s ! coefficient of dissipation. If not specified give just null().
    real(kind = 8), optional, intent(in) :: alpha_opt_v ! coefficient of dissipation. If not specified give just null().
    
    real(kind = 8) :: cosangle, sinangle
    real(kind = 8) :: cartmat2cartel3(3, 3) ! rotation matrix from the material coordinate system to the elemental coordinate system.
    real(kind = 8) :: cartmat2cartel6(6, 6) ! tranformation matrix for the strain in Voight form from the material coordinate system to the elemental coordinate system.

    this%cscriptmat = orthotropicmaterial(cmatv) ! here the orthotropic material type is adopted.

    this%density = density

    this%thickness = thickness
    
    if (present(angledeg_opt)) then

       this%angledeg = angledeg_opt
       
       this%angle = this%angledeg*deg2rad
       
       cosangle = dcos(this%angle)

       sinangle = dsin(this%angle)
    
       this%e1mat = cosangle*e1+sinangle*e2 ! e1mat expresed in the basis of the elemental coordinate system.
                 
       this%e2mat =-sinangle*e1+cosangle*e2 ! e2mat expresed in the basis of the elemental coordinate system.
                 
       this%e3mat = e3 ! e3mat expresed in the basis of the elemental coordinate system.
       
       cartmat2cartel3 = curva2curvb3matrix(this%e1mat, this%e2mat, this%e3mat, e1, e2, e3)
       
       cartmat2cartel6 = curva2curvb6matrix(cartmat2cartel3)
       
       this%cscript = matmul(transpose(cartmat2cartel6), matmul(this%cscriptmat, cartmat2cartel6)) ! transformation of the constitutive law in Voight form into the elemental coordinate system.
       
    else

       this%angledeg = 0.0d0

       this%angle = 0.0d0

       this%e1mat = e1
          
       this%e2mat = e2
                 
       this%e3mat = e3

       this%cscript = this%cscriptmat

    end if
    
    if (present(alpha_opt_v)) then
       this%alpha_v = alpha_opt_v
    else
       this%alpha_v = 0.0d0
    end if      

    if (present(alpha_opt_s)) then
       this%alpha_s = alpha_opt_s
    else
       this%alpha_s = 0.0d0
    end if      
    
    return
    
  end subroutine initializationsinglelayermaterial
    
end module class_single_layer_material
