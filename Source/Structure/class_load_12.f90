module class_load_12

  implicit none

  !> 6 DOF load class
  type :: load_12

     integer :: node               !< node to apply load to
     real(kind = 8) :: spatial(6)  !< spatial load: force (1:3), moment(4:6)
     real(kind = 8) :: material(6) !< material load: force (1:3), moment(4:6)
          
   contains
     
  end type load_12

contains
  
end module class_load_12
