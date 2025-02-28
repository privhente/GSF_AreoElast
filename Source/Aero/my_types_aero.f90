module my_types_aero
  implicit none

  type tangent_matrices3_3x3
     real(kind = 8) :: kvs1(3,3) = 0.0d0, kvs2(3,3) = 0.0d0, kvs3(3,3)=0.0d0
  end type tangent_matrices3_3x3
  
contains

end module my_types_aero