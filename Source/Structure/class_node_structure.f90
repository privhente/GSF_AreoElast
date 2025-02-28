module class_node_structure
  
  use my_constants_structure, only: i_33, e1, e2, e3
  use my_math_structure, only: cross, skew, outer
  
  implicit none
  
  type node_structure
     
     integer :: localid
     integer :: globalid
     
  end type node_structure

end module class_node_structure
