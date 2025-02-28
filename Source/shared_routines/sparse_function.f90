module sparse_function
  
  implicit none
  
    type :: single_sparse_matrix
      integer :: iCallCounter !< current assembly handle
      integer :: nRows, nColumns, nCoords !< number of rows and columns 
    
      integer, allocatable, dimension(:) :: coo_rows, coo_cols !< coordinate format cofficient positions
      integer, allocatable, dimension(:) :: csr_rowIndices, csr_columns !< compressed sparse row format cofficient positions
    
      integer, allocatable, dimension(:) :: handle_row_indices, handle_col_indices !< row and column indices for handle assembly
      integer, allocatable, dimension(:) :: handle_value_indices !< indices to the CSR value vector
    
      integer, allocatable, dimension(:) :: handle_indices !< start indices for handles in handle_???_indices vectors
    
      real(kind = 8), dimension(:), allocatable :: coo_values !< coordinate format cofficient values
      real(kind = 8), dimension(:), allocatable :: csr_values !< compressed sparse row format cofficient values
      
    end type
    
  contains
  
  subroutine sparse_initialize(this,nRows,nColumns)
    implicit none
    class(single_sparse_matrix) :: this
    integer :: nColumns, nRows
    
    this%nColumns = nColumns
    this%nRows = nRows
    
    if (allocated(this%coo_values)) deallocate(this%coo_values)
    if (allocated(this%coo_cols)) deallocate(this%coo_cols)
    if (allocated(this%coo_rows)) deallocate(this%coo_rows)
    if (allocated(this%handle_row_indices)) deallocate(this%handle_row_indices)
    if (allocated(this%handle_col_indices)) deallocate(this%handle_col_indices)
    if (allocated(this%handle_indices)) deallocate(this%handle_indices)
    if (allocated(this%handle_value_indices)) deallocate(this%handle_value_indices)
    if (allocated(this%csr_values)) deallocate(this%csr_values)
    if (allocated(this%csr_rowindices)) deallocate(this%csr_rowindices)
    if (allocated(this%csr_columns)) deallocate(this%csr_columns)

    this%iCallCounter = 0
    this%nCoords = 0
    
  end subroutine  
  
  subroutine create_rowcol_format(this, row_indices, col_indices, matrix)
    implicit none
    class(single_sparse_matrix) :: this
    
    integer, dimension(:), intent(in) :: row_indices, col_indices
    integer :: iOffset, row, nNonZeros, col
    integer :: nRows, nColumns
    real(kind = 8), dimension(:,:), intent(in) :: matrix

    this%iCallCounter =  this%iCallCounter + 1

    nRows = size(row_indices)
    nColumns = size(col_indices)
  
    ! first make room for data
    iOffset = this%nCoords + 1
      
    ! search for zeroes
    nNonZeros = count(matrix .ne. 0.0d0)
      
    ! reallocate for nonzeroes
    call reallocI(this%coo_rows, nNonZeros, this%nCoords)
    call reallocI(this%coo_cols, nNonZeros, this%nCoords)
    call reallocF(this%coo_values, nNonZeros, this%nCoords)
       
    ! reallocate assembly indices
    call reallocI(this%handle_row_indices, nNonZeros, this%nCoords)
    call reallocI(this%handle_col_indices, nNonZeros, this%nCoords)
    
    ! create new index
    call reallocI(this%handle_indices, 1, this%iCallCounter)
    this%handle_indices(this%iCallCounter) = this%nCoords + 1
    
    ! put nonzero values into coordinate format storage
    do row = 1, nRows
      do col = 1, nColumns
        if ( matrix(row, col) /= 0.0d0 ) then
          
          this%coo_rows(iOffset) = row_indices(row)
          this%coo_cols(iOffset) = col_indices(col)
          this%coo_values(iOffset) = matrix(row, col)
          
          this%handle_row_indices(iOffset) = row
          this%handle_col_indices(iOffset) = col
          
          iOffset = iOffset + 1
        endif
      end do
    end do
    
    this%nCoords = this%nCoords + nNonZeros
    return
    
  end subroutine 
  
! =================================================================================================================
subroutine create_csr_format(this)
    implicit none
    class(single_sparse_matrix) :: this    ! finalize handle indices
    
    integer :: iReadOffset, iSearch
    integer :: row, rowIndex, rowIndexDup, col, temp
    integer, allocatable, dimension(:) :: rowValueIndices, csrValueIndices, csrValueIndicesSorted
    integer, allocatable, dimension(:) :: sortingTracker, rowCounts, rowIndicesDupli, columnsDedup
    logical :: bDeduplicated
    
    real(kind = 8), dimension(:), allocatable :: valuesDedup
    real(kind = 8) :: tempf
    
    logical :: bConverged
    
    call reallocI(this%handle_indices, 1, this%iCallCounter)
    this%handle_indices(this%iCallCounter + 1) = this%nCoords + 1

    ! allocate counters for row elements
    allocate(rowCounts(this%nRows))
    
    ! count how often each row occurs
    rowCounts(:) = 0
    
    do iReadOffset = 1, this%nCoords
      row = this%coo_rows(iReadOffset)
      rowCounts(row) = rowCounts(row) + 1
    end do
    
    ! create row indices vector
    allocate(rowIndicesDupli(this%nRows + 1) )
    
    rowIndicesDupli(1) = 1
    do row = 1, this%nRows
      rowIndicesDupli(row + 1) = rowIndicesDupli(row) + rowCounts(row)
    end do
    
    ! count how many deduplicated elements per row
    rowCounts(:) = 0
    
    ! allocate deduplicated columns and values
    allocate(columnsDedup(this%nCoords))
    
    allocate(valuesDedup(this%nCoords))
    allocate(rowValueIndices(this%nCoords))

    ! accumulate elements in same position, remove zeroes and deduplicate
    ! create mapping from handle to csr value vector
    do iReadOffset = 1, this%nCoords
      bDeduplicated = .FALSE.
      
      row = this%coo_rows(iReadOffset)
        
      if(this%coo_values(iReadOffset) /= 0.0d0) then
        
        do iSearch = rowIndicesDupli(row), rowIndicesDupli(row) + rowCounts(row) - 1
          
          if (this%coo_cols(iReadOffset) .eq. columnsDedup(iSearch) )  then
            valuesDedup(iSearch) = valuesDedup(iSearch) + this%coo_values(iReadOffset)
            
            rowValueIndices(iReadOffset) = iSearch
            
            bDeduplicated = .TRUE.
            exit
          end if
        end do
      else
        ! flag this input as irrelevant
        rowValueIndices(iReadOffset) = 1
        bDeduplicated = .TRUE.
      end if
        
      if (.not. bDeduplicated) then
        col = rowIndicesDupli(row) + rowCounts(row)
        
        columnsDedup(col) = this%coo_cols(iReadOffset)
        valuesDedup(col) = this%coo_values(iReadOffset)
        
        rowValueIndices(iReadOffset) = col
        
        rowCounts(row) = rowCounts(row) + 1
      end if    
    end do
      
    ! copy to compressed sparse row format
    if (allocated(this%csr_columns)) then
      deallocate(this%csr_columns)
    end if
    
    if (allocated(this%csr_values)) then
      deallocate(this%csr_values)
    end if
    
    if (allocated(this%csr_rowIndices)) then
      deallocate(this%csr_rowIndices)
    end if

    allocate(this%csr_columns(sum(rowCounts)))
    allocate(this%csr_values(sum(rowCounts)))
    allocate(this%csr_rowIndices(this%nRows + 1))
    allocate(csrValueIndices(this%nCoords))
  
    
    this%csr_rowIndices(1) = 1
    do row = 1, this%nRows
      ! sort columns
      rowIndex = this%csr_rowIndices(row)
      rowIndexDup = rowIndicesDupli(row)
      
      this%csr_columns(rowIndex : rowIndex + rowCounts(row) - 1) = columnsDedup(rowIndexDup : rowIndexDup + rowCounts(row) - 1)
      this%csr_values (rowIndex : rowIndex + rowCounts(row) - 1) = valuesDedup (rowIndexDup : rowIndexDup + rowCounts(row) - 1)
      
      this%csr_rowIndices(row + 1) = this%csr_rowIndices(row) + rowCounts(row)
      
      do iReadOffset = 0, rowCounts(row) - 1
        csrValueIndices(rowIndexDup + iReadOffset) = rowIndex + iReadOffset
      end do
      
    end do
    
    allocate(sortingTracker(size(this%csr_columns)))
    
    do iReadOffset = 1, size(this%csr_columns)
      sortingTracker(iReadOffset) = iReadOffset
    end do
    
    ! sort each row by the columns index in ascending order. keep track of swaps by shuffledPositions
    
    do row = 1, this%nRows
      ! naive reordering: sort as long as there are swaps
      bConverged = .FALSE.
      do while (.not. bConverged)
        bConverged = .TRUE.
        ! check every element for correct order
        ! loop until one before last element (-2)
        do col = this%csr_rowIndices(row), this%csr_rowIndices(row + 1) - 2
          ! if order is not correct, do a swap
          if (this%csr_columns(col) > this%csr_columns(col + 1) ) then
            temp = this%csr_columns(col)
            this%csr_columns(col) = this%csr_columns(col + 1)
            this%csr_columns(col + 1) = temp
            
            tempf = this%csr_values(col)
            this%csr_values(col) = this%csr_values(col + 1)
            this%csr_values(col + 1) = tempf
            
            ! also swap shuffled positions
            temp = sortingTracker(col)
            sortingTracker(col) = sortingTracker(col + 1)
            sortingTracker(col + 1) = temp
            
            bConverged = .FALSE.
          endif
        end do
      end do
    end do
    
    allocate(csrValueIndicesSorted(size(this%csr_columns)))
    do iReadOffset = 1, size(this%csr_columns)
      csrValueIndicesSorted(sortingTracker(iReadOffset)) = iReadOffset
    end do
    
    ! allocate mapping vector to CSR values
    if (allocated(this%handle_value_indices)) then
      deallocate(this%handle_value_indices)
    end if
    
    allocate(this%handle_value_indices(this%nCoords))
    do iReadOffset = 1, this%nCoords
      this%handle_value_indices(iReadOffset) = csrValueIndicesSorted(csrValueIndices(rowValueIndices(iReadOffset)))
    end do
    
    if (allocated(this%coo_rows)) deallocate(this%coo_rows)
    if (allocated(this%coo_cols)) deallocate(this%coo_cols)
    if (allocated(this%coo_values)) deallocate(this%coo_values)
    if (allocated(rowValueIndices)) deallocate(rowValueIndices)
    if (allocated(csrValueIndices)) deallocate(csrValueIndices)
    if (allocated(sortingTracker)) deallocate(sortingTracker)
    if (allocated(csrValueIndicesSorted)) deallocate(csrValueIndicesSorted)
    
    if (allocated(rowCounts)) deallocate(rowCounts)
    if (allocated(rowIndicesDupli)) deallocate(rowIndicesDupli)
    if (allocated(columnsDedup)) deallocate(columnsDedup)
    if (allocated(valuesDedup)) deallocate(valuesDedup)
    
    return
    
  end subroutine 

  !> Reallocate an integer vector, only works for making the array larger
  
  !> @param[in] array the array to enlarge
  !> @param[in] iAddedSize number of elements to add
  !> @param[in] iCurrentSize number of elements currently in use
  subroutine reallocI(array, iAddedSize, iCurrentSize)
    integer, dimension(:), allocatable :: temp
    integer, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: iAddedSize, iCurrentSize
    integer :: sizearray
    
    if (.NOT. ALLOCATED(array)) then
      sizearray = 0
    else
      sizearray = size(array)
    end if
    
    if((sizearray < iCurrentSize + iAddedSize) .OR. .NOT. ALLOCATED(array)) then
      !allocate(temp(iCurrentSize + iAddedSize + 1))
      allocate(temp(sizearray + 16*1024*1024))
      if(sizearray > 0 .AND. ALLOCATED(array)) then
        temp(1:sizearray) = array
      endif
      call move_alloc(temp, array) 
    endif
    return
  end subroutine
    
  !> Reallocate a double precision vector, only works for making the array larger
  
  !> @param[in] array the array to enlarge
  !> @param[in] iAddedSize number of elements to add
  !> @param[in] iCurrentSize number of elements currently in use
  subroutine reallocF(array, iAddedSize, iCurrentSize)
    real(kind = 8), dimension(:), allocatable :: temp
    real(kind = 8), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: iAddedSize, iCurrentSize
    integer :: sizearray
    
    if (.NOT. ALLOCATED(array)) then
      sizearray = 0
    else
      sizearray = size(array)
    end if
  
    if((sizearray < iCurrentSize + iAddedSize) .OR. .NOT. ALLOCATED(array)) then
      allocate(temp(sizearray + 16*1024*1024))
      if(sizearray > 0 .AND. ALLOCATED(array)) then
        temp(1:sizearray) = array
      endif
      call move_alloc(temp, array)
    endif
    return
  end subroutine
  
end module sparse_function
