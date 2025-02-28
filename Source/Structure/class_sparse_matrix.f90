module class_sparse_matrix
  implicit none

  !> @brief management class for sparse matrix data
  !
  !> Enables straightforward assembly of constant matrices and efficient scheme for updating
  !> sparse matrix coefficients for iterative algorithms.
  !>
  !> Block matrices are passed into the class via assemblyIndexed().
  !> Data passed on to these routines is temporarily stored in coordinate format until a call to initSparseStorage() transforms it to compressed sparse row (CSR) format.
  !>
  !> The first assembly determines which matrix positions are nonzero, so use random initialization for non-constant matrices.
  !> After the first assembly, call initSparseStorage() to get CSR data format.
  !>
  !> For the following iteration clear the assembler state by calling startAssembly().
  !> All the subsequent assemblies of the matrix will directly be transferred to the CSR structure.
  !>
  !> WARNING: the order of assembly has to be the same every time. Index positions are tracked by call counts
  type sparse_matrix
    integer, private :: iCallCounter=0 !< current assembly handle
    integer :: nRows=0, nColumns=0, nCoords=0 !< number of rows and columns 
    
    integer, allocatable, private, dimension(:) :: coo_rows, coo_cols !< coordinate format cofficient positions
    integer, allocatable, dimension(:) :: csr_rowIndices, csr_columns !< compressed sparse row format cofficient positions
    
    integer, allocatable, private, dimension(:) :: handle_row_indices, handle_col_indices !< row and column indices for handle assembly
    integer, allocatable, private, dimension(:) :: handle_value_indices !< indices to the CSR value vector
    
    integer, allocatable, private, dimension(:) :: handle_indices !< start indices for handles in handle_???_indices vectors
    
    real(kind = 8), dimension(:), private, allocatable :: coo_values !< coordinate format cofficient values
    real(kind = 8), dimension(:), allocatable :: csr_values !< compressed sparse row format cofficient values
     
  contains
     
    procedure :: initSparseStorage
    procedure :: assemblyIndexed
    procedure :: startAssembly
    procedure :: overWriteValue
    
  end type sparse_matrix
     
  interface sparse_matrix
    module procedure create_sparse_matrix
  end interface
  
  private reallocI, reallocF
  
  contains
  
!> Overwrite values in sparse matrix; Values will only be overwritten, if there exists an entry at the called position in the sparse matrix
  
!> @param[in] row_indices indices for rows
!> @param[in] col_indices indices for columns
!> @param[in] m the block matrix to assemble
subroutine overWriteValue(this, row_indices, col_indices, m)
    implicit none
    class(sparse_matrix) :: this
    integer, dimension(:), intent(in) :: row_indices, col_indices
    integer :: k, i, j, size_vec_range
    integer :: r1, r2
    integer, allocatable, dimension(:) :: vec_range
    real(kind = 8), intent(in) :: m(:,:)
    
  ! Loop over all row indices
    do i = 1,size(row_indices)
      ! getting size of row indicx i; column and values of row index
      size_vec_range = size(this%csr_values(this%csr_RowIndices(row_indices(i)):this%csr_RowIndices(row_indices(i)+1)-1))
      allocate(vec_range(size_vec_range))
      vec_range = 0
      ! get indices of columns of row i
      !vec_range(:) = reshape( (/ this%csr_RowIndices(row_indices(i)):this%csr_RowIndices(row_indices(i)+1)-1 /),(/size_vec_range/))
      
      r1 = this%csr_RowIndices(row_indices(i))
      r2 = this%csr_RowIndices(row_indices(i)+1)-1
      vec_range(:) = (/r1,r2/)

      ! loop over all columns to find same index and replace inside csr_values
      do k = 1,size_vec_range
        do j = 1,size(col_indices)
          if (this%csr_Columns(vec_range(k)) == col_indices(j)) then
            this%csr_values(vec_range(k)) = m(i, j)
            exit
          end if
        end do
      end do
      deallocate(vec_range)
    end do
    
end subroutine overWriteValue

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
    
  !> Constructor for sparse matrix class.
  !
  !> @param[in] nRows number of rows
  !> @param[in] nColumns number of columns
  type(sparse_matrix) function create_sparse_matrix(nRows, nColumns )
    integer, intent(in) :: nRows, nColumns   
  
    create_sparse_matrix%nColumns = nColumns
    create_sparse_matrix%nRows = nRows
    create_sparse_matrix%iCallCounter = 0
    create_sparse_matrix%nCoords = 0
    return
  end function

  !> Clear matrix data and assembler state.
  subroutine startAssembly(this)
    class(sparse_matrix) :: this
    this%csr_values = 0.0d0
    this%iCallCounter = 0
  end subroutine
  
  !> @brief Direct assembly for constant sparse matrices
  !
  !> Directly assemble block matrix into the matrix. This uses temporary storage in coordinate format.
  !> Also takes care of removing zero entries to conserve memory space.
  !
  !> @param[in] row_indices indices for rows
  !> @param[in] col_indices indices for columns
  !> @param[in] m the block matrix to assemble
  !> @param[in] bAllvalues if true, all values of m will be written in sparse matrix, i.e. rows and columns of zero values, too.
  subroutine assemblyIndexed(this, row_indices, col_indices, m, bAllvalues)
    implicit none
    class(sparse_matrix) :: this
    integer, dimension(:), intent(in) :: row_indices, col_indices
      
    real(kind = 8), dimension(:,:), intent(in) :: m
      
    integer :: iOffset, row, nNonZeros, col, assembly_index
    integer :: nRows, nColumns
    
    logical, optional, intent(in) :: bAllvalues
    logical :: bAllvals
    
    bAllvals = .FALSE.
    if (present(bAllvalues)) then
      if (bAllvalues) then
        bAllvals = .True.
      endif
    endif
    
    this%iCallCounter =  this%iCallCounter + 1
    
    if (allocated(this%handle_value_indices)) then
      ! if this field is allocated, we already ran deduplication and the sparse data field is ready to go
      ! first find out where our indices start
      do assembly_index = this%handle_indices(this%iCallCounter), this%handle_indices(this%iCallCounter + 1) - 1
          this%csr_values(this%handle_value_indices(assembly_index)) = this%csr_values(this%handle_value_indices(assembly_index)) + m(this%handle_row_indices(assembly_index), this%handle_col_indices(assembly_index))
      end do
      ! break routine
      return
    endif
    
    nRows = size(row_indices)
    nColumns = size(col_indices)
  
    ! first make room for data
    iOffset = this%nCoords + 1
      
    ! search for zeroes
    if (bAllvals) then
      nNonZeros = size(m,1)*size(m,2)
    else
      nNonZeros = count(m /= 0.0d0)
    end if
    
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
        
        if (bAllvals) then
            this%coo_rows(iOffset) = row_indices(row)
            this%coo_cols(iOffset) = col_indices(col)
            this%coo_values(iOffset) = m(row, col)
            this%handle_row_indices(iOffset) = row
            this%handle_col_indices(iOffset) = col
            iOffset = iOffset + 1
        else
          if ( m(row, col) /= 0.0d0 ) then
            this%coo_rows(iOffset) = row_indices(row)
            this%coo_cols(iOffset) = col_indices(col)
            this%coo_values(iOffset) = m(row, col)
            this%handle_row_indices(iOffset) = row
            this%handle_col_indices(iOffset) = col
            iOffset = iOffset + 1
          endif
        endif
        
      end do
    end do
    
    this%nCoords = this%nCoords + nNonZeros
    return
  end subroutine

  !> @brief Initialize compressed sparse row format sparse storage.
  !
  !> Accumulates coordinates and coefficients previously recorded to the appropriate places into a compressed sparse row structure.
  !> This is done by first deduplicating the coordinate format data and accumulating coefficients
  !>
  !> Finally the compressed sparse row structure is assembled.
  subroutine initSparseStorage(this)
    implicit none
    class(sparse_matrix) :: this

    integer :: iReadOffset, iSearch
      
    integer :: row, rowIndex, rowIndexDup, col, temp
    logical :: bDeduplicated
    integer, allocatable, dimension(:) :: rowValueIndices, csrValueIndices, csrValueIndicesSorted, sortingTracker, rowCounts, rowIndicesDupli, columnsDedup
    real(kind = 8), dimension(:), allocatable :: valuesDedup
    
    real(kind = 8) :: tempf
    
    logical :: bConverged
    
    ! finalize handle indices
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
      
    ! accumulate elements in same position, remove zeros and deduplicate
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
    allocate(this%csr_columns(sum(rowCounts)) )
    allocate(this%csr_values(size(this%csr_columns) ) )
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

    allocate(this%handle_value_indices(this%nCoords))
    
    do iReadOffset = 1, this%nCoords
      this%handle_value_indices(iReadOffset) = csrValueIndicesSorted(csrValueIndices(rowValueIndices(iReadOffset)))
    end do
    
    deallocate(this%coo_rows)
    deallocate(this%coo_cols)
    deallocate(this%coo_values)
    deallocate(rowValueIndices)
    deallocate(csrValueIndices)
    deallocate(sortingTracker)
    deallocate(csrValueIndicesSorted)
    
    deallocate(rowCounts)
    deallocate(rowIndicesDupli)
    deallocate(columnsDedup)
    deallocate(valuesDedup)
      
    return
  end subroutine
    
end module class_sparse_matrix
