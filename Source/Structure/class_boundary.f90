module class_boundary

  implicit none

  !> Boundary conditions for constraints with 12 DOF
  type :: boundary_12
     integer :: constraint12                      !< constraint12 index to apply BC to
   contains
     
   end type boundary_12
   
  !> Boundary conditions for constraints with 6 DOF
  type :: boundary_6
     integer :: constraint6           !< constraint6 index to apply BC to
   contains
     
  end type boundary_6

contains

end module class_boundary
