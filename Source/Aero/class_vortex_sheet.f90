module class_vortex_sheet

  use class_node3_aero
  use class_vortex_segment
  use class_vortex_ring
  use sparse_function, only: single_sparse_matrix, create_rowcol_format, sparse_initialize
  implicit none
  
  type :: vortex_sheet

     integer :: nnodes
     integer :: nrings
     integer :: nsegments
     integer :: ncoordinates
     integer :: nnodes1
     integer :: nnodes2
     integer, allocatable :: indicesq(:), indicesc(:), indicesv(:), indicesp(:)
     real(kind = 8) :: cutoff
     real(kind = 8) :: diffusion_coefficient = 0.0d0
     type(node3_aero), allocatable :: nodes(:)
     type(vortex_ring), allocatable :: rings(:)
     type(vortex_segment), allocatable :: segments(:)
     real(kind = 8), allocatable :: ringscirculations(:), oldringscirculations(:)
     real(kind = 8), allocatable :: segmentscirculations(:), oldsegmentscirculations(:)
     real(kind = 8), allocatable :: q_t(:), v_t(:)
     type(single_sparse_matrix) :: sparse_Trs

   contains
     procedure :: rings2segments_matrix
     
  end type vortex_sheet
  contains

subroutine rings2segments_matrix(this)  
    implicit none
    class(vortex_sheet) :: this
    integer :: i, arr_i(1), arr_j(1)
    integer :: convert_job_coo_csr(8), info
    real(kind = 8) :: arr_val(1,1)
    
    convert_job_coo_csr    = 0 ! the matrix in the CSR format is converted to the coordinate format;
    convert_job_coo_csr(1) = 2 ! the matrix in the coordinate format is converted to the CSR format, and the column indices in CSR representation are sorted in the increasing order within each row.
    convert_job_coo_csr(2) = 1 ! one-based indexing for the matrix in CSR format is used.
    convert_job_coo_csr(3) = 1 ! one-based indexing for the matrix in coordinate format is used.
    convert_job_coo_csr(6) = 0 ! all arrays acsr, ja, ia are filled in for the output storage.
    
    ! initialize sparse matrix
    call sparse_initialize(this%sparse_Trs,this%nsegments,this%nrings)
    
    !< constructing the sparse matrix in COO format
    do i = 1,this%nsegments
        arr_i   = i
        if (this%segments(i)%adjacency(1) /= 0) then
            arr_j   = this%segments(i)%adjacency(1)
            arr_val = 1.0d0
            call create_rowcol_format(this%sparse_Trs, arr_i, arr_j, arr_val)
        end if
        if (this%segments(i)%adjacency(2) /= 0) then
            arr_j   = this%segments(i)%adjacency(2)
            arr_val = -1.0d0
            call create_rowcol_format(this%sparse_Trs, arr_i, arr_j, arr_val)
        end if
    end do
    
    !< storing the sparse matrix in CSR format
    allocate(this%sparse_Trs%csr_values(this%sparse_Trs%ncoords))
    allocate(this%sparse_Trs%csr_columns(this%sparse_Trs%ncoords))
    allocate(this%sparse_Trs%csr_rowIndices(this%nsegments+1))
    
    ! convert coo format to csr format
    !call mkl_scsrcoo(job, n, acsr, ja, ia, nnz, acoo, rowind, colind, info)
    convert_job_coo_csr(5) = this%sparse_Trs%ncoords ! maximum number of the non-zero elements allowed if job(1)=0.
    call mkl_dcsrcoo(convert_job_coo_csr, this%nsegments, this%sparse_Trs%csr_values, &
                  this%sparse_Trs%csr_columns, this%sparse_Trs%csr_rowIndices, this%sparse_Trs%ncoords, &
                  this%sparse_Trs%coo_values(1:this%sparse_Trs%ncoords), this%sparse_Trs%coo_rows(1:this%sparse_Trs%ncoords), this%sparse_Trs%coo_cols(1:this%sparse_Trs%ncoords), info)
    deallocate(this%sparse_Trs%coo_values)
    deallocate(this%sparse_Trs%coo_rows)
    deallocate(this%sparse_Trs%coo_cols)
    
    return
end subroutine rings2segments_matrix

end module class_vortex_sheet
