module class_unbounded_vortex_sheet
  use class_vortex_sheet
  use my_constants_aero,  only : eps_val
  implicit none
  
  type, extends(vortex_sheet) :: unbounded_vortex_sheet
     !type(vortex_segment), allocatable :: wake_edge_segments(:)
     integer, allocatable :: enodes(:)
     integer, allocatable :: erings(:)
     integer :: surface ! to which surface it belongs
     integer :: nnodest ! t from time
     integer :: nringst
     integer :: nnodess ! s from space
     integer :: nringss
     integer :: nproperty
     integer :: nsteps
     integer :: nstepstotal, nrowstotal, nrowstoconsider, tcounter
     integer, allocatable :: segments_ID_first_row(:), rings_ID_first_row(:), rings_ID_first_row_bounded_sheet(:)
     integer, allocatable :: segments_t(:)
     real(kind = 8), allocatable :: q_t_first_row(:), temp_ringscirculations(:)
     logical :: bool_CutWake = .FALSE. , bool_reallocate = .TRUE.
  contains
    procedure :: computeconnectivity
    procedure :: velocity => unboundedvortexsheetvelocity
    procedure :: setproperty => unboundedvortexsheetsetproperty
    procedure :: setgeometry => unboundedvortexsheetsetgeometry
    procedure :: reallocate_arrays
  end type unbounded_vortex_sheet

  contains

 !< Subroutine to reallocate wake position and circulation vector
  subroutine reallocate_arrays(this)
    implicit none
    class(unbounded_vortex_sheet) :: this
    real(kind = 8), allocatable :: temp_array(:)
    
    !< position vector
    allocate(temp_array(3*this%nnodess*(this%tcounter+1)))
    temp_array(:) = this%q_t(1:3*this%nnodess*(this%tcounter+1))
          
    deallocate(this%q_t)
    allocate(this%q_t(3*this%nnodess*(this%tcounter+1)))
    this%q_t(:) = temp_array
     
    !< ring circulation vector - in case this has to be done like above, consider alse to reduce transformation matrix for segments!
    this%ringscirculations(this%tcounter*this%nringss+1:this%nrings) = 0.0d0
    
    deallocate(temp_array)
    this%bool_reallocate = .FALSE.
    
    return    
  end subroutine reallocate_arrays
  
!< Subroutine to compute wake-node connectivity
  subroutine computeconnectivity(this)
    implicit none
    class(unbounded_vortex_sheet) :: this
    integer :: i, j, k, l
    
    do i = 1, this%nringst
       do j = 1, this%nringss
          k = (i-1)*this%nringss+j
          l = (i-1)*this%nnodess+j
          this%rings(k)%connectivity(:) = [ l+this%nnodess+1, l+this%nnodess, l, l+1]
       end do
    end do
    
    return
  end subroutine computeconnectivity

!< Subroutine to set and actualize wake geometry
  subroutine unboundedvortexsheetsetgeometry(this, boolean_first_Wake_row_only_opt)
    implicit none
    
    class(unbounded_vortex_sheet), intent(inout) :: this
    integer :: i, nsegments
    logical :: boolean_first_Wake_row_only
    logical, optional :: boolean_first_Wake_row_only_opt
    real(kind = 8) :: temp_q_t_rings(12), temp_q_t_segments(6)

    boolean_first_Wake_row_only = .False.
    if (present(boolean_first_Wake_row_only_opt)) boolean_first_Wake_row_only = boolean_first_Wake_row_only_opt
    
!< Loop over all segments to assign new segment geometry and segment circulations
    nsegments = 0
    if (boolean_first_Wake_row_only .eqv. .False.) then
      if (allocated(this%segments_t)) nsegments = size(this%segments_t)
      !$OMP PARALLEL DO & 
      !$OMP SHARED(nsegments) &
      !$OMP PRIVATE(i)   
      do i = 1, nsegments
        temp_q_t_segments=this%q_t(this%segments(this%segments_t(i))%indicesq)
        !call this%segments(this%segments_t(i))%setgeometry(this%q_t(this%segments(this%segments_t(i))%indicesq))
        call this%segments(this%segments_t(i))%setgeometry(temp_q_t_segments)
      end do
      !$OMP END PARALLEL DO  
    else ! if boolean_first_Wake_row_only = .True. update only the first row segments but not segment bounding to 2nd wake row!
      nsegments = size(this%segments_ID_first_row)
      do i = 1, nsegments
        temp_q_t_segments=this%q_t(this%segments(this%segments_ID_first_row(i))%indicesq)
        !call this%segments(this%segments_ID_first_row(i))%setgeometry(this%q_t(this%segments(this%segments_ID_first_row(i))%indicesq))
        call this%segments(this%segments_ID_first_row(i))%setgeometry(temp_q_t_segments)
      end do
    end if
    
    !< update ring geometry (coordinates, plane area, -normal, -tangent), surface velocities and free-field velocity at ring control point
    !< but only for the first wake row
    !!$OMP PARALLEL DO & 
    !!$OMP PRIVATE(i)   
    do i = 1, this%nringss
      !< update geometry and surface velocity
      temp_q_t_rings=this%q_t(this%rings(i)%indicesq)
      !call this%rings(i)%setgeometry(this%q_t(this%rings(i)%indicesq))
      call this%rings(i)%setgeometry(temp_q_t_rings)
    end do
    !!$OMP END PARALLEL DO  
   
    return
  end subroutine unboundedvortexsheetsetgeometry
  
!< Subroutine to set and actualize wake properties (geometry and circulation)
  subroutine unboundedvortexsheetsetproperty(this)
    implicit none
    
    class(unbounded_vortex_sheet), intent(inout) :: this
    integer :: i, nsegments
    real(kind = 8), allocatable :: circulations(:)
    integer, allocatable :: indices(:), choices(:)
    real(kind = 8) :: temp_q_t_segments(6)
    
!< sparse multiplication to calculate segment circulations: Gamma = T*G   
    allocate(circulations(this%nrings))
    circulations = this%ringscirculations
    call mkl_dcsrmv('N', this%nsegments, this%nrings, 1.0d0, 'G  F  ', &
      this%sparse_Trs%csr_values, this%sparse_Trs%csr_columns, &
      this%sparse_Trs%csr_rowIndices(1:size(this%sparse_Trs%csr_rowIndices) - 1), &
      this%sparse_Trs%csr_rowIndices(2:size(this%sparse_Trs%csr_rowIndices)), &
      circulations, 0.0d0, this%segmentscirculations)

!< find only segments with non-zero circulation and only consider them
    nsegments = count(abs(this%segmentscirculations) .gt. eps_val )
    allocate(indices(this%nsegments))
    allocate(choices(nsegments))
    indices = merge( 0, [ ( i, i = 1, this%nsegments ) ], abs(this%segmentscirculations) .lt. eps_val )
    if (allocated(this%segments_t)) deallocate(this%segments_t)

    allocate(this%segments_t(nsegments))
    this%segments_t = pack( indices, indices .ne. 0 )
  
!< Loop over all segments to assign new segment geometry and segment circulations
    !$OMP PARALLEL DO & 
    !$OMP SHARED(nsegments) &
    !$OMP PRIVATE(i)   
    do i = 1, nsegments
      temp_q_t_segments=this%q_t(this%segments(this%segments_t(i))%indicesq)
      !call this%segments(this%segments_t(i))%setgeometry(this%q_t(this%segments(this%segments_t(i))%indicesq))
      call this%segments(this%segments_t(i))%setgeometry(temp_q_t_segments)
      !< segment circulation
      this%segments(this%segments_t(i))%circulation = this%segmentscirculations(this%segments_t(i))
    end do
    !$OMP END PARALLEL DO  
    
    return
  end subroutine unboundedvortexsheetsetproperty
  
!< Function to compute velocity due to wake segments 
  function unboundedvortexsheetvelocity(this, targetpoint, cutoff, time) result(v)
    implicit none
    
    class(unbounded_vortex_sheet), intent(inout) :: this
    real(kind = 8), intent(in) :: targetpoint(3)
    real(kind = 8), intent(in) :: cutoff, time
    real(kind = 8) :: circulation
    integer :: i, size_wake_segments
    real(kind = 8) :: v(3)
    
    v(:) = 0.0d0
    
    size_wake_segments = 0
    if (allocated(this%segments_t)) size_wake_segments = size(this%segments_t)

    !!$OMP PARALLEL DO & 
    !!$OMP SHARED(size_wake_segments,cutoff,targetpoint) &
    !!$OMP PRIVATE(i)   
    do i = 1, size_wake_segments
      circulation = this%segments(this%segments_t(i))%circulation
      this%segments(this%segments_t(i))%time   = time
      this%segments(this%segments_t(i))%cutoff = sqrt ( cutoff**2 + 4.0d0*1.25643d0*this%diffusion_coefficient*abs(circulation)*time )
      v = v + this%segments(this%segments_t(i))%velocity(targetpoint)
    end do
    !!$OMP END PARALLEL DO
    
    return
  end function unboundedvortexsheetvelocity
    
end module class_unbounded_vortex_sheet
