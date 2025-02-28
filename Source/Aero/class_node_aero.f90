module class_node_aero
  
  use my_constants_aero, only: i_33, e1, e2, e3
  use my_math_aero, only: cross, skew, outer
  
  implicit none
  
  type node_aero
     
     integer :: localid
     integer :: globalid
     
  end type node_aero

end module class_node_aero
