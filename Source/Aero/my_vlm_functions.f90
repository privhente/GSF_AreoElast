module my_vlm_functions
  use my_constants_aero, only: pi, i_33
  use my_math_aero, only: cross, skew, outer
  implicit none
  
  contains

!< Function to calculate the vorticity induced velocity using the exact solution of velocity for vortex segments with two points, see Predictman page 46
  function segment_velocity(point1, point2, circulation, cutoff, targetpoint)
    implicit none
    
    real(kind = 8) :: segment_velocity(3)
    real(kind = 8), dimension(3), intent(in) :: point1, point2, targetpoint
    real(kind = 8), intent(in) :: circulation, cutoff
    real(kind = 8), dimension(3) :: r1, r2
    real(kind = 8) :: nr1, nr2, nl
    
    r1 = targetpoint-point1
    r2 = targetpoint-point2
    nl = norm2(r2-r1)
    nr1 = norm2(r1)
    nr2 = norm2(r2)
    
    !< analytical solution of vorticity induced velocity
    segment_velocity = circulation*(nr1+nr2)*cross(r1, r2)/(4.0d0*pi*(nr1*nr2*(nr1*nr2+dot_product(r1, r2)) + (cutoff*nl)**2))
       
    return
  end function segment_velocity
  
!< Subroutine for writing character line from file to an array/vector
  subroutine read_line_to_array(file_unit, vec, var_type, ireason)
    implicit none
    
    integer, intent(in):: file_unit
    integer, intent(inout):: ireason
    integer :: k, var_type
    
    character(len = 1) :: strline
    character(len=:), allocatable :: temp_strLine
    
    real(kind = 8) :: val
    real(kind = 8), allocatable, intent(inout) :: vec(:)
    
    logical :: boolean_valfound = .FALSE.

    temp_strLine = ''
    if (allocated(vec)) deallocate(vec)
    k = 0
    do
      read(file_unit,'(a)',advance='NO',iostat=ireason) strline
      
      if (k .eq. 0) then
        if (ireason .ne. 0) then
            ireason = -1
            exit
        end if
      end if
      
      if (strline .ne. ' ') then
        boolean_valfound = .TRUE.
        temp_strLine = trim(adjustl(temp_strLine)) // trim(adjustl(strline))
      else
        if (boolean_valfound) then
          if (var_type == 0) then
            read(temp_strLine, '(f30.15)' )  val
            call append_real_vec(vec, val) 
          end if
          
        end if
        temp_strLine = ''
        boolean_valfound = .FALSE.
      end if
      
      if (IS_IOSTAT_END(ireason)) then
        ireason = 0
        exit
      end if

      k = k + 1
    end do
          
    return
  end subroutine read_line_to_array

!< subroutine for appending real values to an array
  subroutine append_real_vec(vec, val)
    implicit none
    real(kind = 8), allocatable, intent(inout) :: vec(:)
    real(kind = 8), intent(in) :: val
    real(kind = 8), allocatable :: temp_vec(:)
    integer :: num_add
    
    if (.not. allocated(vec)) then
      allocate(vec(1))
      num_add = 0
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
  
!< Subroutine to compute and store segments, connectivities, segment and ring adjancies
  subroutine find_segments(nnodes1, nnodes2, segmentsconnectivities, segmentsadjacencies, ringsadjacencies)
    implicit none

    integer, intent(in) :: nnodes1, nnodes2
    integer, intent(inout) :: segmentsconnectivities(:,:)
    integer, intent(inout) :: segmentsadjacencies(:,:)
    integer, intent(inout) :: ringsadjacencies(:,:)
    integer, allocatable, dimension(:) :: ones, seqa, seqb, seq1, seq2, seq3, seq4
    integer :: nsegments, nelements
    integer :: li, ui
    integer :: i, j, k, l, insegments
      
    nsegments = nnodes1*(nnodes2-1)+(nnodes1-1)*nnodes2
    nelements = (nnodes1-1)*(nnodes2-1)
  
    allocate(ones(nnodes1-1))
    allocate(seqa(nnodes1-1))
    allocate(seqb(nnodes1-1))
    allocate(seq1(nnodes1-1))
    allocate(seq2(nnodes1-1))
    allocate(seq3(nnodes1-1))
    allocate(seq4(nnodes1-1))

    !< finding conectivities for the segments
    insegments = 1
    do i = 1, nnodes1-1
       do j = 1, nnodes2
          k = nnodes1*(j-1)+i
          segmentsconnectivities(insegments, :) = [k, k+1]
          insegments = insegments+1
       end do
    end do 

    do i = 1, nnodes2-1
       do j = 1, nnodes1
          k = nnodes1*(i-1)+j
          segmentsconnectivities(insegments, :) = [k, k+nnodes1] 
          insegments = insegments+1
       end do
    end do

  !< finding adjacencies for the segments
  !< first edge
    do i = 1, nnodes1-1
       k = nnodes2*(i-1)+1
       segmentsadjacencies(k, :) = [i, 0] 
    end do

  !< intern between first and third edges
    do i = 2, nnodes2-1
       do j = 1, nnodes1-1
          k = nnodes2*(j-1)+i
          l = (nnodes1-1)*(i-2)+j
          segmentsadjacencies(k, :) = [l+nnodes1-1, l] 
       end do
    end do
 
  !< third edge
    do i = 1, nnodes1-1
       k = nnodes2*(i-1)+nnodes2
       segmentsadjacencies(k, :) = [0, (nnodes1-1)*(nnodes2-2)+i] 
    end do

  !< second edge
    do i = 1, nnodes2-1
       k = (nnodes1-1)*nnodes2+nnodes1*i
       segmentsadjacencies(k, :) = [(nnodes1-1)*i, 0] 
    end do

  !< intern between second and fourth edges
    do i = 2, nnodes1-1
       do j = 1, nnodes2-1
          k = (nnodes1-1)*nnodes2+nnodes1*(j-1)+i
          l = (nnodes1-1)*(j-1)+i-1
          segmentsadjacencies(k, :) = [l, l+1] 
       end do
    end do
  
  !< fourth edge
    do i = 1, nnodes2-1
       k = (nnodes1-1)*nnodes2+nnodes1*(i-1)+1
       segmentsadjacencies(k, :) = [0, (nnodes1-1)*(i-1)+1] 
    end do

  !< finding the segments that belong to each ring
    ones(:) = 1
    do i = 1, nnodes1-1
       seqa(i) = (i-1)*nnodes2+1
       seqb(i) = (nnodes1-1)*nnodes2+i
    end do

    do i = 1, nnodes2-1
       li = (i-1)*(nnodes1-1)+1
       ui = i*(nnodes1-1)
     
       seq1 = seqa+ones* (i-1)*1
       seq2 = seqb+ones*((i-1)*nnodes1+1)
       seq3 = seqa+ones*((i-1)*1+1)
       seq4 = seqb+ones* (i-1)*nnodes1
     
       ringsadjacencies(li:ui, 1) = seq1
       ringsadjacencies(li:ui, 2) = seq2
       ringsadjacencies(li:ui, 3) = seq3
       ringsadjacencies(li:ui, 4) = seq4
    end do

  !< deallocate unused variables
    deallocate(ones)
    deallocate(seqa)
    deallocate(seqb)
    deallocate(seq1)
    deallocate(seq2)
    deallocate(seq3)
    deallocate(seq4)
  
    return
  end subroutine find_segments

!< Function to calculate the vorticity induced velocity using the exact solution of velocity for vortex segments with two points, see Predictman page 46
  function tangent_segment_velocity(point1, point2, circulation, cutoff, targetpoint)
    implicit none
    
    real(kind = 8), dimension(3), intent(in) :: point1, point2, targetpoint
    real(kind = 8), intent(in) :: cutoff
    real(kind = 8) :: tangent_segment_velocity(3,6)
    real(kind = 8) :: nr1, nr2, nl, r1(3), r2(3), hv(3), dv, dhv(3,6), ddv(6)
    real(kind = 8) :: circulation
    
    r1 = targetpoint-point1
    r2 = targetpoint-point2
    nl = norm2(r2-r1)
    nr1 = norm2(r1)
    nr2 = norm2(r2)
    
    hv = (nr1+nr2)*cross(r1,r2)
    dv = (nr1*nr2*(nr1*nr2+dot_product(r1, r2)) + (cutoff*nl)**2)
    
    dhv(1:3,1:3) = outer(cross(r1,r2),r1)/nr1 - (nr1+nr2)*skew(r2)
    dhv(1:3,4:6) = outer(cross(r1,r2),r2)/nr2 + (nr1+nr2)*skew(r1)
    
    ddv(1:3) = 2*(nr2**2)*r1 + nr1*nr2*r2 + dot_product(r1,r2)*nr2*r1/nr1 - 2.0d0*(cutoff**2)*(r2-r1)
    ddv(4:6) = 2*(nr1**2)*r2 + nr1*nr2*r1 + dot_product(r1,r2)*nr1*r2/nr2 + 2.0d0*(cutoff**2)*(r2-r1)

    !< analytical solution of vortex induced velocity with unit circulation
    tangent_segment_velocity = circulation/(4.0d0*pi) * (dv*dhv - outer(hv,ddv)) / (dv**2)
    
    return
  end function tangent_segment_velocity  
end module my_vlm_functions
  