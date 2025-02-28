module my_fsi_functions
  
  implicit none
  contains
  
  subroutine random_array(array)
    USE ieee_arithmetic  ! available in  module load gcc/9.
    implicit none
    real(kind = 8), allocatable, dimension(:), intent(inout) :: array
    real(kind = 8) :: r
    r = ieee_value(r,  ieee_positive_inf)
    array = r
  end subroutine random_array
  
!< find position of values in an array
  function findpositioninarray(array_s,array_val) result(inz)
    implicit none
    integer, intent(in) :: array_s(:), array_val(:)
    integer, allocatable :: inz(:)
    integer :: i, j
    
    allocate(inz(size(array_val)))
    inz = 0
    do i = 1,size(array_val)
      do j = 1,size(array_s)  
        if (array_val(i) .eq. array_s(j)) then
          inz(i) = j
          exit
        end if
      end do
    end do
    return
  end function findpositioninarray
  
!< function for appending values to an integer  
  subroutine append_int_vec(vec, val, flag)
    implicit none
    integer, allocatable, intent(inout) :: vec(:)
    integer, intent(in) :: val, flag
    integer, allocatable :: temp_vec(:)
    integer :: num_add
    logical :: boolean_statement
    
    if (.not. allocated(vec)) then
      allocate(vec(1))
      vec(1) = -1
      num_add = 0
      vec = 0
    else
      num_add = 1
    end if
    
    boolean_statement = .TRUE. ! append values
    if (flag .eq. 1) boolean_statement = all(vec .ne. val) ! append values only, if not exist already
    
    if (boolean_statement) then
      allocate(temp_vec(size(vec) + num_add))
      temp_vec(1:size(vec)) = vec
      temp_vec(size(vec) + num_add) = val
      deallocate(vec)
      allocate(vec(size(temp_vec)))
      vec = temp_vec
    end if
    
    return
  end subroutine append_int_vec

!< function for appending values to an integer  
  subroutine append_real_vec(vec, val)
    implicit none
    real(kind = 8), allocatable, intent(inout) :: vec(:)
    real(kind = 8), intent(in) :: val
    real(kind = 8), allocatable :: temp_vec(:)
    integer :: num_add
    
    if (.not. allocated(vec)) then
      allocate(vec(1))
      num_add = 0
      vec = 0
    else
      num_add = 1
    end if
    
    allocate(temp_vec(size(vec) + num_add))
    temp_vec(1:size(vec)) = vec
    temp_vec(size(vec) + num_add) = val
    deallocate(vec)
    allocate(vec(size(temp_vec)))
    vec = temp_vec
    
    return
  end subroutine append_real_vec

!< function for appending values to an integer  
  subroutine append_real_vec3(vec, vec3)
    implicit none
    real(kind = 8), allocatable, intent(inout) :: vec(:,:)
    real(kind = 8), intent(in) :: vec3(3)
    real(kind = 8), allocatable :: temp_vec(:,:)
    integer :: num_add
    
    if (.not. allocated(vec)) then
      allocate(vec(1,3))
      num_add = 0
      vec = 0
    else
      num_add = 1
    end if
    
    allocate(temp_vec(size(vec,1) + num_add,3))
    temp_vec(1:size(vec,1),1:3) = vec(1:size(vec,1),1:3)
    temp_vec(size(vec,1) + num_add,1:3) = vec3(1:3)
    deallocate(vec)
    allocate(vec(size(temp_vec,1),3))
    vec = temp_vec
    
    return
  end subroutine append_real_vec3
  
!< function for weighting function
  function weight_factor(L,R,int_flag)
    implicit none
    real(kind = 8) :: L,R, epsilon = 1.0d-6
    real(kind = 8) :: weight_factor
    integer :: int_flag

    !< weighting function
    if (int_flag .eq. 1) then ! bulb function
      weight_factor = exp(-1.0d0/(1.0d0-(L/R)**2))/exp(-1.0d0)
    else if (int_flag .eq. 2) then ! Gaussian function
      weight_factor = exp(log(epsilon)*(L/R)**2)
    end if
    
    return
  end function weight_factor
  
  !< subroutine to determine the determinant of a 3x3 matrix 
  function det(A)
    implicit none
    real(kind = 8), dimension(3,3) :: A
    real(kind = 8) :: det

    det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(3,1)*A(2,2)*A(1,3) - A(3,2)*A(2,3)*A(1,1) - A(3,3)*A(2,1)*A(1,2) 
    
    return
  end function det
  
end module my_fsi_functions