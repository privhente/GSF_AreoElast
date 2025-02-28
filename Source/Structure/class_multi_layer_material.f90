module class_multi_layer_material

  use class_single_layer_material
  
  implicit none

  type :: multi_layer_material

     integer :: nplys ! number of plys.
     integer :: plysequence ! name of the plysequence
     type(single_layer_material), allocatable :: plys(:)
     real(kind = 8) :: thickness ! total thickness of the multi-layer material.
     real(kind = 8) :: masterjacobean
     real(kind = 8), allocatable :: mastermasterjacobean(:)  ! jacobean to perform the layer-wise integration.
     real(kind = 8), allocatable :: boundaries(:) ! boundaries of the plys in natural coordinates.
     real(kind = 8), allocatable :: masterboundaries(:) ! boundaries of the plys in natural coordinates.
     
   contains

     procedure :: initialization => initializationmultilayermaterial
     
!!$     procedure :: boundaryfinder  
!!$     procedure :: inwhichlayer
          
  end type multi_layer_material
  
contains

  subroutine initializationmultilayermaterial(this)

    implicit none
    
    class(multi_layer_material), intent(inout) :: this
    
    integer :: i

    this%thickness = 0.0d0
    
    do i = 1, this%nplys

       this%thickness = this%thickness+this%plys(i)%thickness
       
    end do

    this%masterjacobean = 0.5d0*this%thickness
    
    this%boundaries(1) = -0.5d0*this%thickness

    this%masterboundaries(1) = -1.0d0
    
    do i = 1, this%nplys

       this%mastermasterjacobean(i) =                          this%plys(i)%thickness      /this%thickness
       
       this%boundaries(i+1)         = this%boundaries(i)      +this%plys(i)%thickness
       
       this%masterboundaries(i+1)   = this%masterboundaries(i)+this%plys(i)%thickness*2.0d0/this%thickness
       
    end do    

    return
    
  end subroutine initializationmultilayermaterial
      
end module class_multi_layer_material
