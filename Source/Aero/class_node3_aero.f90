module class_node3_aero
  
  use class_node_aero
  
  implicit none

  type, extends(node_aero) :: node3_aero
     
     integer :: coordinates(3)
     real(kind = 8) :: q_0(3)
     real(kind = 8) :: q_t(3), v_t(3)
     integer :: indicesq(3), indicesv(3)
     real(kind = 8) :: vinfinity(3) = 0.0d0
     
     ! new variables for fsi
     real(kind = 8), allocatable :: Tn(:,:) ! transformation matrix for aero nodes to structural nodes (for FSI)
     integer, allocatable :: nodes_st(:)    ! node list of structural nodes for aero
     integer, allocatable :: indicesq_st_global(:), indicesv_st_global(:)
   contains
     
  end type node3_aero
  
contains

end module class_node3_aero
