module my_math_aero
  
  implicit none
  
contains
  
  function outer(a, b)

    implicit none

    real(kind = 8), intent(in) :: a(:), b(:)

    real(kind = 8) :: outer(size(a), size(b))
    
    integer :: i, j
    
    forall(i = 1:size(a), j = 1:size(b))

       outer(i, j) = a(i)*b(j)

    end forall
    
    return
    
  end function outer

  function cross(a, b)
    
    implicit none

    real(kind = 8), dimension(3) :: cross
    real(kind = 8), dimension(3), intent(in) :: a, b
    
    cross(1) = a(2)*b(3)-a(3)*b(2)
    cross(2) = a(3)*b(1)-a(1)*b(3)
    cross(3) = a(1)*b(2)-a(2)*b(1)
    
    return
    
  end function cross

  function skew(a)

    implicit none

    real(kind = 8), dimension(3, 3) :: skew
    real(kind = 8), dimension(3), intent(in) :: a

    skew(:, :) = 0
    
    skew(1, 2) =-a(3)
    skew(1, 3) = a(2)
    skew(2, 3) =-a(1)

    skew(2, 1) = a(3)
    skew(3, 1) =-a(2)
    skew(3, 2) = a(1)
 
    return

  end function skew

  function diag(a)

    implicit none

    real(kind = 8), intent(in) :: a(:)
    real(kind = 8) :: diag(size(a), size(a))

    integer :: i
    
    diag(:, :) = 0.0d0

    do i = 1, size(a)
    
       diag(i, i) = a(i)

    end do

    return

  end function diag
  
  function eye(n) result(i_nn)

    implicit none
    
    integer, intent(in) :: n
    real(kind=8) :: i_nn(n, n)  
    
    integer :: i  
    
    i_nn(:,:) = 0.0d0
    
    do i = 1, n
       
       i_nn(i, i) = 1.0d0

    end do
    
    return
    
  end function eye

  function vec6mat9symm(vec) result(mat)

    implicit none
    
    real(kind = 8) :: vec(6), mat(3, 3)
    
    mat(1, 1) = vec(1)
    mat(2, 2) = vec(2)
    mat(3, 3) = vec(3)
    
    mat(2, 3) = vec(4)
    mat(1, 3) = vec(5)
    mat(1, 2) = vec(6)
    
    mat(3, 2) = vec(4)
    mat(3, 1) = vec(5)
    mat(2, 1) = vec(6)
    
    return
    
  end function vec6mat9symm

  function vec10mat16symm(vec) result(mat)

    implicit none
    
    real(kind = 8) :: vec(10), mat(4, 4)
    
    mat(1, 1) = vec( 1)
    mat(2, 2) = vec( 2)
    mat(3, 3) = vec( 3)
    mat(4, 4) = vec( 4)
    
    mat(3, 4) = vec( 5)
    mat(2, 4) = vec( 6)
    mat(1, 4) = vec( 7)
    mat(1, 3) = vec( 8)
    mat(1, 2) = vec( 9)
    mat(2, 3) = vec(10)

    mat(4, 3) = vec( 5)
    mat(4, 2) = vec( 6)
    mat(4, 1) = vec( 7)
    mat(3, 1) = vec( 8)
    mat(2, 1) = vec( 9)
    mat(3, 2) = vec(10)
           
    return
    
  end function vec10mat16symm
  
  function disdzfz(x, y, fx, fy, dzfz, tolerance_arg)

    implicit none
    
    real(kind = 8), intent(in) :: x(:), y(:), fx, fy, dzfz(:)
    real(kind = 8), intent(in), optional :: tolerance_arg
    
    real(kind = 8) :: disdzfz(size(x))
    
    real(kind = 8) :: delta(size(x)), sqndelta

    real(kind = 8) :: tolerance = 1.0d-10

    if(present(tolerance_arg)) tolerance = tolerance_arg
    
    delta = y-x
    
    sqndelta = dot_product(delta, delta)
    
    if (dsqrt(sqndelta) < tolerance) then
       
       disdzfz = dzfz
       
    else
       
       disdzfz = dzfz+((fy-fx)-dot_product(dzfz, delta))*delta/sqndelta
       
    end if
    
    return

  end function disdzfz

  function dydisdzfz(x, y, fx, fy, dzfz, dyfy, dzdzfz, tolerance_arg)

    implicit none
   
    real(kind = 8), intent(in) :: x(:), y(:), fx, fy, dzfz(:), dyfy(:), dzdzfz(:, :)
    real(kind = 8), intent(in), optional :: tolerance_arg
    
    real(kind = 8)  :: dydisdzfz(size(x), size(x))

    real(kind = 8) :: delta(size(x)), sqndelta, i_nn(size(x), size(x))

    real(kind = 8) :: tolerance = 1.0d-10

    if(present(tolerance_arg)) tolerance = tolerance_arg
    
    delta = y-x
    
    sqndelta = dot_product(delta, delta)     
    
    if (dsqrt(sqndelta) < tolerance) then
       
       dydisdzfz = 0.5d0*dzdzfz
       
    else
       
       i_nn = eye(size(x))
       
       dydisdzfz = 0.5d0*dzdzfz &
                  +(outer(delta, dyfy-dzfz)-0.5d0*matmul(outer(delta, delta), dzdzfz))/sqndelta &
                  +(fy-fx-dot_product(dzfz, delta))*(i_nn/sqndelta-2.0d0*outer(delta, delta)/sqndelta**2)
       
    end if
    
    return
    
  end function dydisdzfz
  
end module my_math_aero
