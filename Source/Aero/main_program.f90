program main_program
  use class_model_aero, only: model_aero
  implicit none
  type(model_aero) :: the_model_aero
  call main_aero(the_model_aero)
end program main_program