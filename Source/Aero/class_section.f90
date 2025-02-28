module class_section

  implicit none

  type :: section

     integer :: nsegments
     integer :: nrings
     integer, allocatable :: segments(:) 
     integer, allocatable :: rings(:)
!     real, allocatable :: circulations(:)
     
  end type section
  
end module class_section
