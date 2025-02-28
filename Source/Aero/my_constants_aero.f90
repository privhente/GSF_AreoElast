module my_constants_aero
  implicit none

  real(kind = 8), parameter :: pi = 3.141592653589793d0
  real(kind = 8), parameter :: o_33(3, 3) = reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/),(/3, 3/)) 
  real(kind = 8), parameter :: i_33(3, 3) = reshape((/1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/),(/3, 3/))
  real(kind = 8), parameter :: e1(3) = [1.0d0, 0.0d0, 0.0d0]
  real(kind = 8), parameter :: e2(3) = [0.0d0, 1.0d0, 0.0d0]
  real(kind = 8), parameter :: e3(3) = [0.0d0, 0.0d0, 1.0d0]
  real(kind = 8), parameter :: eps_val = 1.0d-16
  
end module my_constants_aero
