module class_nullspace
  
  implicit none
  
  type :: cnode
    integer, allocatable :: indices_temp_q(:)
  end type cnode
  
  type :: clist
    type(cnode), allocatable :: node(:)
    integer :: rank
    integer, allocatable :: constraints(:)
    integer, allocatable :: indicesn(:)
    integer, allocatable :: nodes6(:), nodes12(:), nodes(:)
    integer, allocatable :: H(:,:)
    integer :: n_nodes12, n_nodes6
  end type clist
   
  type :: nsi
      integer :: rank
      integer, allocatable :: indicesn(:), indicesq(:)
  end type nsi

  type :: nullspace
     type(nsi), allocatable :: ns_6(:), ns_12(:), ns_add6(:) 
     type(clist), allocatable :: list(:)
     integer, allocatable :: list12(:,:), list6(:,:), addlist6(:)
     integer, allocatable :: constraint_list(:,:)
  end type nullspace
  

end module class_nullspace
