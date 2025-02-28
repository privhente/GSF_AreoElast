module class_load_6

  implicit none
  
  !> 6 DOF load class
  type :: load_6

     integer :: node              !< node to apply load to
     real(kind = 8) :: spatial(3) !< spatial load
     real(kind = 8) :: material   !< material load
          
   contains
     
  end type load_6

contains
  
end module class_load_6
