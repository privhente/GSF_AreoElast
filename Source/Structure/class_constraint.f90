module class_constraint

  use my_constants_structure, only: i_33

  implicit none
  
  !> base class for constraints
  type :: constraint
     
     integer :: nodes(2)                !< node indices for this constraint, may be one or two
     real(kind = 8):: phi(3, 2) = 0.0d0 , dir(3) = 0.0d0, dir_gl(3) = 0.0d0 !< one shift vector per node, one vector for direction
     
     integer, allocatable :: indicesg(:), indicesq(:), indicesn(:)
     
     integer :: rankcount    !< rank of this constraint
     integer :: coordinates  !< number of nodes involved in this constaint
     logical :: stiffness    !< indicates whether this constraint has a stiffness term
     logical :: boolean_q    !< indicates whether this constraint is depending on q and/or s
     logical :: boolean_v    !< indicates whether this constraint is depending on q and/or s
     real(kind = 8) :: amplitude_translation, amplitude_translation_vel, amplitude_rotation, amplitude_rotation_vel
     real(kind = 8) :: deltat
     
  end type constraint
  
end module class_constraint
