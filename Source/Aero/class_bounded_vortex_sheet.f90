module class_bounded_vortex_sheet

  use class_vortex_sheet
  use class_section

  implicit none
  
  type, extends(vortex_sheet) :: bounded_vortex_sheet
     integer :: nsections !< number of sections of the bounded vortex sheet (no wake sheet)
     type(section), allocatable :: sections(:)
     real(kind = 8), allocatable :: DGamma(:,:), DG(:,:), DGamma_dense(:,:)
     integer, allocatable :: indices_qv_aero(:)
     logical :: boolean_lifting_surface = .FALSE.
    contains
      procedure :: setkinematic => boundedvortexsheetsetkinematic
      procedure :: setproperty => boundedvortexsheetsetproperty
      procedure :: velocity => boundedvortexsheetvelocity
    end type bounded_vortex_sheet

  contains

!< Subroutine to set and update bounded vortex sheet (ring and segment) geometry
  subroutine boundedvortexsheetsetkinematic(this)
    implicit none
    
    class(bounded_vortex_sheet), intent(inout)  :: this  
    integer :: i
    real(kind = 8) :: temp_q_t_segments(6),temp_q_t_rings(12), temp_v_t_rings(12)

!!$OMP PARALLEL DO & 
!!$OMP PRIVATE(i)
    do i = 1, this%nnodes
       this%nodes(i)%q_t(:) = this%q_t(this%nodes(i)%indicesq)
    end do
!!$OMP END PARALLEL DO

    !< update ring geometry (coordinates, plane area, -normal, -tangent), surface velocities and free-field velocity at ring control point
!!$OMP PARALLEL DO & 
!!$OMP PRIVATE(i)
    do i = 1, this%nrings
      !< update geometry and surface velocity
      temp_q_t_rings=this%q_t(this%rings(i)%indicesq)
      call this%rings(i)%setgeometry(temp_q_t_rings)
      !call this%rings(i)%setgeometry(this%q_t(this%rings(i)%indicesq))
      
      temp_v_t_rings=this%v_t(this%rings(i)%indicesv)
      !call this%rings(i)%setsurfvelocity(this%v_t(this%rings(i)%indicesv))
      call this%rings(i)%setsurfvelocity(temp_v_t_rings)
    end do
!!$OMP END PARALLEL DO
    
    !< update segment geometry
!!$OMP PARALLEL DO & 
!!$OMP PRIVATE(i)
    do i = 1, this%nsegments
      temp_q_t_segments=this%q_t(this%segments(i)%indicesq)
      !call this%segments(i)%setgeometry(this%q_t(this%segments(i)%indicesq)) ! segment indices q in global segment direction
      call this%segments(i)%setgeometry(temp_q_t_segments) ! segment indices q in global segment direction
    end do
!!$OMP END PARALLEL DO    
    return
  end subroutine boundedvortexsheetsetkinematic

!< Subroutine to set and update bounded vortex sheet (ring and segment) circulation  
  subroutine boundedvortexsheetsetproperty(this)
    implicit none
    
    class(bounded_vortex_sheet), intent(inout) :: this
    integer :: i
        
!!$OMP PARALLEL DO & 
!!$OMP PRIVATE(i)
    do i = 1, this%nrings
       this%rings(i)%oldcirculation = this%oldringscirculations(i) ! ring circulation from previous steps
       this%rings(i)%circulation    = this%ringscirculations(i)    ! update ring circulations
    end do
!!$OMP END PARALLEL DO
    
    !< transformation of ring circulations to segment circulations
    call mkl_dcsrmv('N', this%nsegments, this%nrings, 1.0d0, 'G  F  ', &
      this%sparse_Trs%csr_values, this%sparse_Trs%csr_columns, &
      this%sparse_Trs%csr_rowIndices(1:size(this%sparse_Trs%csr_rowIndices) - 1), &
      this%sparse_Trs%csr_rowIndices(2:size(this%sparse_Trs%csr_rowIndices)), &
      this%ringscirculations, 0.0d0, this%segmentscirculations)

    call mkl_dcsrmv('N', this%nsegments, this%nrings, 1.0d0, 'G  F  ', &
      this%sparse_Trs%csr_values, this%sparse_Trs%csr_columns, &
      this%sparse_Trs%csr_rowIndices(1:size(this%sparse_Trs%csr_rowIndices) - 1), &
      this%sparse_Trs%csr_rowIndices(2:size(this%sparse_Trs%csr_rowIndices)), &
      this%oldringscirculations, 0.0d0, this%oldsegmentscirculations)
    
!!$OMP PARALLEL DO & 
!!$OMP PRIVATE(i)
    do i = 1, this%nsegments
      this%segments(i)%oldcirculation = this%oldsegmentscirculations(i)
      this%segments(i)%circulation    = this%segmentscirculations(i)
    end do
!!$OMP END PARALLEL DO    
    return
  end subroutine boundedvortexsheetsetproperty

!< Function to calculate velocity on target point due to bounded vortex sheets
  function boundedvortexsheetvelocity(this, targetpoint, cutoff) result(v)
    implicit none
    
    class(bounded_vortex_sheet), intent(inout) :: this
    real(kind = 8), intent(in) :: targetpoint(3)
    real(kind = 8) :: v(3), cutoff
    integer :: i
    
    v(:) = 0.0d0
    do i = 1, this%nsegments
      this%segments(i)%cutoff = cutoff
      v = v+this%segments(i)%velocity(targetpoint)
    end do

    return
  end function boundedvortexsheetvelocity
  
end module class_bounded_vortex_sheet
