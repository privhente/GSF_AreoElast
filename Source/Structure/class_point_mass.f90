module class_point_mass

  implicit none

  type :: point_mass

     real(kind = 8) :: mass
     integer :: node
     
  end type point_mass
  
end module class_point_mass
