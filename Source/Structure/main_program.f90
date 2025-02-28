! DeSiO main program for structural analysis
program main_program
 use class_model_structure, only: model_structure
  implicit none
  type(model_structure) :: the_model_structure
  call main_structure(the_model_structure)
end program main_program