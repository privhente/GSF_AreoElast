module my_math_structure
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

  function kron(i,j)
    implicit none
    integer :: kron
    integer, intent(in) :: i,j
    
    if (i .eq. j) then
        kron = 1
    else
        kron = 0
    end if
    
    return
  end function kron
  
  function levi(i,j,k)
    implicit none
    integer :: levi
    integer, intent(in) :: i,j,k
    
    
    if (i .eq. 1 .AND. j .eq. 2 .AND. k .eq. 3) then
        levi = 1
    elseif (i .eq. 2 .AND. j .eq. 3 .AND. k .eq. 1) then
        levi = 1
    elseif (i .eq. 3 .AND. j .eq. 1 .AND. k .eq. 2) then
        levi = 1
    elseif (i .eq. 1 .AND. j .eq. 3 .AND. k .eq. 2) then
        levi = -1
    elseif (i .eq. 2 .AND. j .eq. 1 .AND. k .eq. 3) then
        levi = -1
    elseif (i .eq. 3 .AND. j .eq. 2 .AND. k .eq. 1) then
        levi = -1
    else 
        levi = 0
    end if 
 
    return
  end function levi
    
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

  function vec21mat36symm(vec) result(mat)
    implicit none
    real(kind = 8) :: vec(21), mat(6, 6)
    ! diagonal
    mat(1, 1) = vec(1)
    mat(2, 2) = vec(2)
    mat(3, 3) = vec(3)
    mat(4, 4) = vec(4) 
    mat(5, 5) = vec(5)
    mat(6, 6) = vec(6)
    ! upper diagonal
    mat(5, 6) = vec(7)
    mat(4, 6) = vec(8)
    mat(3, 6) = vec(9)
    mat(2, 6) = vec(10)
    mat(1, 6) = vec(11)
    mat(1, 5) = vec(12)
    mat(1, 4) = vec(13)
    mat(1, 3) = vec(14)
    mat(1, 2) = vec(15)
    mat(2, 3) = vec(16)
    mat(3, 4) = vec(17)
    mat(4, 5) = vec(18)
    mat(3, 5) = vec(19)
    mat(2, 5) = vec(20)
    mat(2, 4) = vec(21)
    
    ! lower diagonal
    mat(6, 5) = vec(7)
    mat(6, 4) = vec(8)
    mat(6, 3) = vec(9)
    mat(6, 2) = vec(10)
    mat(6, 1) = vec(11)
    mat(5, 1) = vec(12)
    mat(4, 1) = vec(13)
    mat(3, 1) = vec(14)
    mat(2, 1) = vec(15)
    mat(3, 2) = vec(16)
    mat(4, 3) = vec(17)
    mat(5, 4) = vec(18)
    mat(5, 3) = vec(19)
    mat(5, 2) = vec(20)
    mat(4, 2) = vec(21)
    return
  end function vec21mat36symm

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

  function curva2curvb3matrix(asub1, asub2, asub3, bsup1, bsup2, bsup3) result(curva2curvb3)
!!$ this function computes the transformation from the curvilinear system "a" to the curvilinear system "b".
    implicit none
    real(kind = 8), dimension(3), intent(in) :: asub1, asub2, asub3, bsup1, bsup2, bsup3
!!$ asub1, asub2 and asub3 are the basis of "a".
!!$ bsup1, bsup2 and bsup3 are the co-basis "b".
    real(kind = 8) :: curva2curvb3(3, 3)
!!$ curva2curvb3 is the transformation matrix which tranforms vectors in "a" to vectors in "b".
    curva2curvb3(1, 1) = dot_product(bsup1, asub1)
    curva2curvb3(1, 2) = dot_product(bsup1, asub2)
    curva2curvb3(1, 3) = dot_product(bsup1, asub3)
    curva2curvb3(2, 1) = dot_product(bsup2, asub1)
    curva2curvb3(2, 2) = dot_product(bsup2, asub2)
    curva2curvb3(2, 3) = dot_product(bsup2, asub3)
    curva2curvb3(3, 1) = dot_product(bsup3, asub1)
    curva2curvb3(3, 2) = dot_product(bsup3, asub2)
    curva2curvb3(3, 3) = dot_product(bsup3, asub3)        
    return
  end function curva2curvb3matrix
  
  function curva2curvb6matrix(curva2curvb3) result(curva2curvb6)
!!$ this function computes the transformation for the strain in Voigth form from the curvilinear system "a" to the curvilinear system "b".
    implicit none
    real(kind = 8), intent(in) :: curva2curvb3(3, 3)
!!$ curva2curvb3 is the transformation matrix which tranforms vectors in "a" to vectors in "b".
    real(kind = 8) :: curva2curvb6(6, 6)
!!$ curva2curvb6 is the transformation matrix which tranforms the strain in Voigth form from "a" to "b".    
    curva2curvb6(1, 1) =       curva2curvb3(1, 1)*curva2curvb3(1, 1)
    curva2curvb6(1, 2) =       curva2curvb3(2, 1)*curva2curvb3(2, 1)
    curva2curvb6(1, 3) =       curva2curvb3(3, 1)*curva2curvb3(3, 1)
    curva2curvb6(1, 4) =       curva2curvb3(2, 1)*curva2curvb3(3, 1)
    curva2curvb6(1, 5) =       curva2curvb3(1, 1)*curva2curvb3(3, 1)
    curva2curvb6(1, 6) =       curva2curvb3(1, 1)*curva2curvb3(2, 1)
    curva2curvb6(2, 1) =       curva2curvb3(1, 2)*curva2curvb3(1, 2)
    curva2curvb6(2, 2) =       curva2curvb3(2, 2)*curva2curvb3(2, 2)
    curva2curvb6(2, 3) =       curva2curvb3(3, 2)*curva2curvb3(3, 2)
    curva2curvb6(2, 4) =       curva2curvb3(2, 2)*curva2curvb3(3, 2)
    curva2curvb6(2, 5) =       curva2curvb3(1, 2)*curva2curvb3(3, 2)
    curva2curvb6(2, 6) =       curva2curvb3(1, 2)*curva2curvb3(2, 2)
    curva2curvb6(3, 1) =       curva2curvb3(1, 3)*curva2curvb3(1, 3)
    curva2curvb6(3, 2) =       curva2curvb3(2, 3)*curva2curvb3(2, 3)
    curva2curvb6(3, 3) =       curva2curvb3(3, 3)*curva2curvb3(3, 3)
    curva2curvb6(3, 4) =       curva2curvb3(2, 3)*curva2curvb3(3, 3)
    curva2curvb6(3, 5) =       curva2curvb3(1, 3)*curva2curvb3(3, 3)
    curva2curvb6(3, 6) =       curva2curvb3(1, 3)*curva2curvb3(2, 3)
    curva2curvb6(4, 1) = 2.0d0*curva2curvb3(1, 2)*curva2curvb3(1, 3)
    curva2curvb6(4, 2) = 2.0d0*curva2curvb3(2, 2)*curva2curvb3(2, 3)
    curva2curvb6(4, 3) = 2.0d0*curva2curvb3(3, 2)*curva2curvb3(3, 3)
    curva2curvb6(4, 4) =       curva2curvb3(2, 2)*curva2curvb3(3, 3)+curva2curvb3(3, 2)*curva2curvb3(2, 3)
    curva2curvb6(4, 5) =       curva2curvb3(1, 2)*curva2curvb3(3, 3)+curva2curvb3(3, 2)*curva2curvb3(1, 3)
    curva2curvb6(4, 6) =       curva2curvb3(1, 2)*curva2curvb3(2, 3)+curva2curvb3(2, 2)*curva2curvb3(1, 3)
    curva2curvb6(5, 1) = 2.0d0*curva2curvb3(1, 1)*curva2curvb3(1, 3)
    curva2curvb6(5, 2) = 2.0d0*curva2curvb3(2, 1)*curva2curvb3(2, 3)
    curva2curvb6(5, 3) = 2.0d0*curva2curvb3(3, 1)*curva2curvb3(3, 3)
    curva2curvb6(5, 4) =       curva2curvb3(2, 1)*curva2curvb3(3, 3)+curva2curvb3(3, 1)*curva2curvb3(2, 3)
    curva2curvb6(5, 5) =       curva2curvb3(1, 1)*curva2curvb3(3, 3)+curva2curvb3(3, 1)*curva2curvb3(1, 3)
    curva2curvb6(5, 6) =       curva2curvb3(1, 1)*curva2curvb3(2, 3)+curva2curvb3(2, 1)*curva2curvb3(1, 3)
    curva2curvb6(6, 1) = 2.0d0*curva2curvb3(1, 1)*curva2curvb3(1, 2)
    curva2curvb6(6, 2) = 2.0d0*curva2curvb3(2, 1)*curva2curvb3(2, 2)
    curva2curvb6(6, 3) = 2.0d0*curva2curvb3(3, 1)*curva2curvb3(3, 2)
    curva2curvb6(6, 4) =       curva2curvb3(2, 1)*curva2curvb3(3, 2)+curva2curvb3(3, 1)*curva2curvb3(2, 2)
    curva2curvb6(6, 5) =       curva2curvb3(1, 1)*curva2curvb3(3, 2)+curva2curvb3(3, 1)*curva2curvb3(1, 2)
    curva2curvb6(6, 6) =       curva2curvb3(1, 1)*curva2curvb3(2, 2)+curva2curvb3(2, 1)*curva2curvb3(1, 2)
    return
  end function curva2curvb6matrix

  end module my_math_structure