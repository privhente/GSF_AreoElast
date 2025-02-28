module solver_functions

  implicit none

  contains

  !> @brief Assemble the iteration matrix from the simulation model
  !
  !> Fill an array with random numbers [1...2]
  !
  !> @param[in] array the array to fill
  subroutine random_array(array)
    USE ieee_arithmetic  ! available in  module load gcc/9.
    implicit none
    real(kind = 8), allocatable, dimension(:), intent(inout) :: array
    real(kind = 8) :: r
    r = ieee_value(r,  ieee_positive_inf)
    array = r
  end subroutine

! subroutine for initialization of global arrays and parameters
  subroutine solver_initialization(this)
    use class_model_structure, only: model_structure
    use solver_variables
    use mkl_dss
    use class_sparse_matrix

    implicit none
    type(model_structure), intent(inout) :: this

    logical :: boolean_abort = .FALSE.
    integer :: ilast

    tnnodes = this%tnnodes
    tncoordinates = this%tncoordinates
    trconstraints = this%trconstraints
    nMatSize = 2*tncoordinates+trconstraints

    allocate(q0 (tncoordinates))
    allocate(q1 (tncoordinates))
    allocate(v0 (tncoordinates))
    allocate(v1 (tncoordinates))
    allocate(p1 (tncoordinates))
    allocate(lambda(trconstraints))
    allocate(lambda0(trconstraints))
    allocate(gravity_vector(tncoordinates))

    sparse_smatrix = sparse_matrix(nMatSize, nMatSize)

    allocate(ml_t(6*tnnodes))
    allocate(sl_t(6*tnnodes))

    ml_t(:) = 0.0d0
    sl_t(:) = 0.0d0

    allocate(indicesq(tncoordinates))
    allocate(indicesv(tncoordinates))
    allocate(indicesg(trconstraints))
    allocate(indicesq_mod(this%tncoordinates12 + this%tncoordinates6))
    allocate(indicesg_stat(trconstraints))

! indicesq - indicesq = body 1,2,... node12_q, body 1 - node6_q, strain8, body 2 - node6_q, strain8, body 3 ..., ...
    do i = 1, tncoordinates
       indicesq(i) = i
       indicesv(i) = tncoordinates+i
    end do

 ! indicesq_mod - for buckling and modal analysis - to avoid strain dofs
    do i = 1, this%nnodes12
      indicesq_mod(12*(i-1)+1:12*(i-1)+12) = this%nodes12(i)%coordinates
    end do
    ilast = this%nnodes12*12
    do i = 1, this%nnodes6
      indicesq_mod(ilast + 6*(i-1)+1 : ilast + 6*(i-1)+6) = this%nodes6(i)%coordinates
    end do

 ! indicesg - constraints
    do i = 1, trconstraints
      indicesg(i) = 2*tncoordinates+i
      indicesg_stat(i) = tncoordinates+i
    end do

! File Format
3232 format(1000000f20.5)
3233 format(1000000i10)

! fill inputs with random data
    print *, "finding non-zero entries in the iteration matrix"

    CALL RANDOM_SEED

    call random_array(lambda)
    call random_array(q0)
    call random_array(q1)
    call random_array(v0)
    call random_array(v1)
    call random_array(ml_t)
    call random_array(sl_t)

! extract used fields from problem matrix
    this%deltat = 1.0d0
    call this%modelboundaries(0.0d0)

! Internal forces, tangent stiffness matrix
    call this%modelfk(q0, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q1, v0, v1)

! computing g and dg for the whole model at q0 and q1.
    call this%modelgdg(q1, q0, sparse_smatrix, indicesg, indicesq, lambda)

    print *, "Initializing sparse iteration matrix..."

    call sparse_smatrix%initSparseStorage()
    print *, "sparse non-zeroes:", size(sparse_smatrix%csr_values)

    ALLOCATE(perm(nMatSize))
    ALLOCATE(perm1(nMatSize-tncoordinates))
    ALLOCATE(perm2(nMatSize-tncoordinates))

  4001 format(i10)

    if(boolean_abort) then
  ! write sparse matrix entries
      open(unit = 4001, file = 'smatrix_rowIndex.txt', action = 'write', status = 'replace')
      write(4001, 4001) sparse_smatrix%csr_rowIndices
      close(4001)

      open(unit = 4002, file = 'smatrix_columns.txt', action = 'write', status = 'replace')
      write(4002, 4001) sparse_smatrix%csr_columns
      close(4002)

    endif

#if B05_ADDON > 1
    call b05_analyze_system(0, 0, sparse_smatrix%csr_values, sparse_smatrix%csr_columns, sparse_smatrix%csr_rowIndices, size(sparse_smatrix%csr_values), size(sparse_smatrix%csr_columns), size(sparse_smatrix%csr_rowIndices))
#endif

    return
  end subroutine solver_initialization

  ! B05: logging routines for the new linear solver, which is implemented in C++.
  subroutine log_res(r_abs, r_rel)
    real(kind = 8), intent(in) :: r_abs, r_rel
    print '(8X,A24,1P2E10.2E2)', 'Squared residual abs/rel', r_abs, r_rel
  end subroutine log_res
  subroutine warn_dense()
    print *, 'WARNING: sparse solve failed! Falling back to dense LAPACK solver'
  end subroutine warn_dense

! Subroutine for the solution procedure using sparse or dense matrices
#if !B05_ADDON
  subroutine solver(nDim, handle, csr_values, csr_columns, csr_rowIndices, rhs_unbalance_force_solver, nRHS, delta_solver, boolean_sparse, statOut, boolean_determinant)
    use class_sparse_matrix
    use solver_variables
    use mkl_dss

    implicit none

    type(MKL_DSS_HANDLE), intent(inout) :: handle

    integer, intent(in) :: nDim, nRHS
    integer, dimension(:), intent(in) :: csr_columns, csr_rowIndices
    integer :: convert_info
    integer, dimension(:), allocatable :: ipiv

    real(kind = 8), intent(in), dimension(:) :: csr_values
    real(kind = 8), intent(in), dimension(:) :: rhs_unbalance_force_solver(:)
    real(kind = 8), intent(inout), dimension(:) :: delta_solver(:)

    real(kind = 8), dimension(:,:), allocatable :: dense_matrix
    logical, intent(in) :: boolean_sparse, boolean_determinant

    integer, parameter :: bufLen = 20
    character*15 :: statIn
    integer :: buff(bufLen) = 0

    real(kind = 8), allocatable, intent(inout) :: statOUt(:)

    if (allocated(rhs_unbalance_force_sol)) deallocate(rhs_unbalance_force_sol)
    allocate(rhs_unbalance_force_sol(nDim))
    rhs_unbalance_force_sol = 0.0d0
    delta_solver            = 0.0d0

    if (allocated(statOut)) then
      deallocate(statOut)
      allocate(statOut(5))
    else
      allocate(statOut(5))
    end if
    statOut = 0.0d0

! solution of newton's iteration
    if (boolean_sparse) then
      error = dss_factor_real(handle, MKL_DSS_POSITIVE_DEFINITE, csr_values)
      error = dss_solve_real(handle, MKL_DSS_DEFAULTS, -rhs_unbalance_force_solver, nRHS, delta_solver)

! Print Out the determinant of the matrix (no statistics for a diagonal matrix)
      if (boolean_determinant) then
        statIn = 'determinant'
        call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
        error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
      end if

! perform solution check: first multiply solution with problem matrix
      call mkl_dcsrmv('N', nDim, nDim, 1.0d0, 'G  F  ', csr_values, csr_columns, csr_rowIndices(1:size(csr_rowIndices) - 1), csr_rowIndices(2:size(csr_rowIndices)), delta_solver, 0.0d0, rhs_unbalance_force_sol)

! Next check the dot product of the deviation
      rhs_unbalance_force_deviation = dsqrt(dot_product(rhs_unbalance_force_sol + rhs_unbalance_force_solver , rhs_unbalance_force_sol + rhs_unbalance_force_solver))
      rhs_unbalance_force_amplitude = dsqrt(dot_product(rhs_unbalance_force_solver, rhs_unbalance_force_solver))
      ! B05: for comparison with new solver. Disable if undesired.
      !call log_res(rhs_unbalance_force_deviation, &
      !     rhs_unbalance_force_deviation / rhs_unbalance_force_amplitude)

      ! check rel. error between residual vector and calculated residual vector
      if (abs(rhs_unbalance_force_deviation / rhs_unbalance_force_amplitude) > 1.0d-6) then
        ! check simple deviation in residual vector
        if (abs(rhs_unbalance_force_deviation) > 1.0d-6) then

          print *, "WARNING: Failed to solve the equations accurately! Falling back to dense solver"
          allocate(dense_matrix(nDim, nDim))
          allocate(ipiv(nDim))
          convert_job(1) = 1
          convert_job(2) = 1
          convert_job(3) = 1
          convert_job(4) = 2

          ! expand sparse matrix
          call mkl_ddnscsr(convert_job, nDim, nDim, dense_matrix, nDim, csr_values, csr_columns, csr_rowIndices, convert_info)

          ! solve dense matrix
          delta_solver = -rhs_unbalance_force_solver
          call dgesv(nDim, 1, dense_matrix, nDim, ipiv, delta_solver, nDim, convert_info)

          deallocate(dense_matrix)
          deallocate(ipiv)
        endif
      end if

    else
        allocate(dense_matrix(nDim, nDim))
        allocate(ipiv(nDim))
        convert_job(1) = 1
        convert_job(2) = 1
        convert_job(3) = 1
        convert_job(4) = 2

        ! expand sparse matrix
        call mkl_ddnscsr(convert_job, nDim, nDim, dense_matrix, nDim, csr_values, csr_columns, csr_rowIndices, convert_info)

        ! solve dense matrix
        delta_solver = -rhs_unbalance_force_solver
        call dgesv(nDim, 1, dense_matrix, nDim, ipiv, delta_solver, nDim, convert_info)

        deallocate(dense_matrix)
        deallocate(ipiv)
    endif

    deallocate(rhs_unbalance_force_sol)
  end subroutine solver
#endif

!< Subroutine to extract sparse value for static case (without velocity dof)
  subroutine extract_csr_Format_From_Sparsematrix(sparsematrix,csr_values, csr_columns,csr_rowIndices, csr_compr_source, indicesRow,indicesCol)
! output indices (rows, columns) and values are reordered new
      use class_sparse_matrix
      implicit none

      type(sparse_matrix), intent(in) :: sparsematrix

      integer, dimension(:), allocatable, intent(inout) :: csr_columns, csr_rowIndices, csr_compr_source
      real(kind = 8), dimension(:), allocatable, intent(inout) :: csr_values
      integer, dimension(:), intent(in) :: indicesRow,indicesCol
      integer :: aa, bb, cc, i, j ,k, r1

      integer, dimension(:), allocatable :: csr_columns_temp, csr_rowIndices_temp, csr_compr_source_temp, columns_n, vec_range
      real(kind = 8), dimension(:), allocatable :: csr_values_temp, values_n

      integer :: size_vec_range
      integer, allocatable, dimension(:) :: csr_Columns_inp, csr_RowIndices_inp
      real(kind = 8), allocatable, dimension(:) :: csr_Values_inp

      allocate(csr_Values_inp(size(sparsematrix%csr_values)))
      allocate(csr_RowIndices_inp(size(sparsematrix%csr_RowIndices)))
      allocate(csr_Columns_inp(size(sparsematrix%csr_Columns)))

      csr_Values_inp     = 0.0d0
      csr_RowIndices_inp = 0
      csr_Columns_inp    = 0

      csr_Values_inp     = sparsematrix%csr_values
      csr_RowIndices_inp = sparsematrix%csr_RowIndices
      csr_Columns_inp    = sparsematrix%csr_Columns

      allocate(csr_Values_temp(size(sparsematrix%csr_values)))
      allocate(csr_Columns_temp(size(sparsematrix%csr_values)))
      allocate(csr_RowIndices_temp(size(sparsematrix%csr_RowIndices)))
      allocate(csr_compr_source_temp(size(sparsematrix%csr_values)))

      csr_Values_temp       = 0.0d0
      csr_Columns_temp      = 0
      csr_RowIndices_temp   = 0
      csr_compr_source_temp = 0

      aa = 0
      cc = 1
      csr_rowIndices_temp(1) = 1
      do i = 1,size(indicesRow)
        size_vec_range = size(csr_Values_inp(csr_RowIndices_inp(indicesRow(i)):csr_RowIndices_inp(indicesRow(i)+1)-1))

        allocate(vec_range(size_vec_range))
        allocate(values_n(size_vec_range))
        allocate(columns_n(size_vec_range))

        vec_range = 0
        values_n  = 0.0d0
        columns_n = 0

        !vec_range = reshape( (/csr_RowIndices_inp(indicesRow(i)):csr_RowIndices_inp(indicesRow(i)+1)-1/),(/size_vec_range/))
        !values_n  = reshape( (/csr_Values_inp(vec_range)/),(/size_vec_range/))
        !columns_n = reshape( (/csr_Columns_inp(vec_range)/),(/size_vec_range/))

        vec_range(1:size_vec_range) = (/ (r1, r1 = csr_RowIndices_inp(indicesRow(i)),csr_RowIndices_inp(indicesRow(i)+1)-1) /)
        values_n(1:size_vec_range)  = csr_Values_inp(vec_range)
        columns_n(1:size_vec_range) = csr_Columns_inp(vec_range)

        bb = 0
        do k = 1,size(columns_n)

          do j = 1,size(indicesCol)
            if (columns_n(k) == indicesCol(j)) then
              aa = aa + 1
              bb = bb + 1
              csr_values_temp(aa)  = values_n(k)
              csr_columns_temp(aa) = j
              csr_compr_source_temp(aa) = vec_range(k)
              exit
            end if
          end do
        end do
        csr_rowIndices_temp(cc+1) = csr_rowIndices_temp(cc) + bb
        cc = cc + 1

        deallocate(vec_range)
        deallocate(values_n)
        deallocate(columns_n)

      end do

      if (allocated(csr_values)) then
        deallocate(csr_values)
      end if
      if (allocated(csr_columns)) then
        deallocate(csr_columns)
      end if
      if (allocated(csr_rowIndices)) then
        deallocate(csr_rowIndices)
      end if
      if (allocated(csr_compr_source)) then
        deallocate(csr_compr_source)
      end if

      allocate(csr_rowIndices(cc))
      allocate(csr_columns(aa))
      allocate(csr_values(aa))
      allocate(csr_compr_source(aa))

      csr_rowIndices   = 0
      csr_columns      = 0
      csr_values       = 0.0d0
      csr_compr_source = 0

      csr_rowIndices(:)   = csr_rowIndices_temp(1:cc)
      csr_columns(:)      = csr_columns_temp(1:aa)
      csr_values(:)       = csr_values_temp(1:aa)
      csr_compr_source(:) = csr_compr_source_temp(1:aa)

      deallocate(csr_Values_temp)
      deallocate(csr_Columns_temp)
      deallocate(csr_RowIndices_temp)
      deallocate(csr_compr_source_temp)
      deallocate(csr_Values_inp)
      deallocate(csr_RowIndices_inp)
      deallocate(csr_Columns_inp)

  end subroutine extract_csr_Format_From_Sparsematrix

  subroutine linSorting_int(inzToSort,arrayToSort)
    implicit none

    logical :: bConverged

    integer :: i, j, k, r1
    integer, allocatable, dimension(:), intent(inout) :: inzToSort
    integer :: inzToCheck, inzToCheck_t

    integer, allocatable, dimension(:), intent(inout) :: arrayToSort
    integer :: valToCheck, valToCheck_t

    allocate(inzToSort(size(arrayToSort)))
    !inzToSort = reshape((/1:size(arrayToSort)/),(/size(arrayToSort)/))
    inzToSort(1:size(arrayToSort)) = (/ (r1, r1=1,size(arrayToSort)) /)

    bConverged = .FALSE.
    i = 1
    k = 0
    do while (i /= size(arrayToSort))
      valToCheck = arrayToSort(i)
      inzToCheck = inzToSort(i)
      do while (k+1 /= size(arrayToSort))
        k = k + 1
        if (arrayToSort(k+1)==valToCheck) then
          k = k + 1
        elseif (arrayToSort(k+1)<valToCheck) then
          valToCheck_t     = arrayToSort(k+1)
          inzToCheck_t     = inzToSort(k+1)
          arrayToSort(k+1) = valToCheck
          inzToSort(k+1)   = inzToCheck
          arrayToSort(i)   = valToCheck_t
          inzToSort(i)     = inzToCheck_t
          k = k - 1
          exit
        end if
        if (k+1 >= size(arrayToSort)) then
          j = 1
          i = i + 1
          k = i - 1
          exit
        end if
        if (k + 1 == size(arrayToSort) .or. (arrayToSort(k+1)==valToCheck)) then
          i  = i + 1
          k  = i - 1
          exit
        end if
      end do
    end do
  end subroutine linSorting_int

  subroutine linSorting_real(inzToSort,arrayToSort)
    implicit none

    logical :: bConverged

    integer :: i, j, k, r1
    integer, allocatable, dimension(:), intent(inout) :: inzToSort
    integer :: inzToCheck, inzToCheck_t

    real(kind = 8), dimension(:), intent(inout) :: arrayToSort
    real(kind =8) :: valToCheck, valToCheck_t

    allocate(inzToSort(size(arrayToSort)))
    !inzToSort = reshape((/1:size(arrayToSort)/),(/size(arrayToSort)/))
    inzToSort(1:size(arrayToSort)) = (/ (r1, r1=1,size(arrayToSort)) /)

    bConverged = .FALSE.
    i = 1
    k = 0
    do while (i /= size(arrayToSort))
      valToCheck = arrayToSort(i)
      inzToCheck = inzToSort(i)
      do while (k+1 /= size(arrayToSort))
        k = k + 1
        if (arrayToSort(k+1)==valToCheck) then
          k = k + 1
        elseif (arrayToSort(k+1)<valToCheck) then
          valToCheck_t     = arrayToSort(k+1)
          inzToCheck_t     = inzToSort(k+1)
          arrayToSort(k+1) = valToCheck
          inzToSort(k+1)   = inzToCheck
          arrayToSort(i)   = valToCheck_t
          inzToSort(i)     = inzToCheck_t
          k = k - 1
          exit
        end if
        if (k+1 >= size(arrayToSort)) then
          j = 1
          i = i + 1
          k = i - 1
          exit
        end if
        if (k + 1 == size(arrayToSort) .or. (arrayToSort(k+1)==valToCheck)) then
          i  = i + 1
          k  = i - 1
          exit
        end if
      end do
    end do
  end subroutine linSorting_real

! Lapack routine to solve generalized Eigenvalue problem
  subroutine get_eigenvalues(eigenvalues,eigenvectors,matrix1,matrix2,nEF,nq,boolean_abort)

    real(kind = 8), intent(in) :: matrix1(:,:), matrix2(:,:)
    real(kind = 8), allocatable, intent(inout) :: eigenvalues(:), eigenvectors(:,:)
    real(kind = 8), allocatable, dimension(:)   :: alpha_real, alpha_imaginary, absz, eigenvalues_temp
    real(kind = 8), allocatable, dimension(:,:) :: vl, vr
    real(kind = 8), allocatable, dimension(:)   :: beta, work
    integer :: info, nq, aa, i, lwork
    integer :: nEF
    integer, allocatable, dimension(:) :: indices_eig, indices_eig_temp, indices_sort_ascent
    logical, intent(inout) :: boolean_abort

    lwork = 8*nq+16

    allocate(work(lwork))
    allocate(alpha_real(nq))
    allocate(alpha_imaginary(nq))
    allocate(vl(nq,nq))
    allocate(vr(nq,nq))
    allocate(beta(nq))

! solving eigenvalue problem with dense matrices: Ax*x*lambda = B*x
    call DGGEV('N','V', nq, matrix1, nq, matrix2, nq, alpha_real, alpha_imaginary, beta, vl, nq, vr, nq, work, lwork, info)

    if (info/=0) then
      print*, 'Eigenvalue analysis failed !!'
      boolean_abort = .True.
      return
    end if

! Sort out imaginary parts of eigenvalues
    allocate(indices_eig_temp(nq))
    allocate(absz(nq))
    aa = 0
    do i = 1,nq
      if (abs(alpha_imaginary(i)) <= 1.0d-8) then
        if (beta(i) /= 0.0) then
          aa = aa + 1
          indices_eig_temp(aa) = i
        end if
      end if
    end do

    allocate(indices_eig(aa))
    allocate(eigenvalues_temp(aa))
    indices_eig(:) = indices_eig_temp(1:aa)
    eigenvalues_temp(:) = (alpha_imaginary(indices_eig) + alpha_real(indices_eig))/beta(indices_eig)

! Sort eigenvalues in ascending order
    call linSorting_real(indices_sort_ascent,eigenvalues_temp)
    do i = 1,size(eigenvalues_temp)
      if (abs(eigenvalues_temp(i)) <= 1.0d-8) then
        eigenvalues_temp(i) = 0.0d0
      end if
    end do
    nEF = aa
    allocate(eigenvectors(nq,nEF))
    allocate(eigenvalues(nEF))
    eigenvectors(:,1:nEF) = vr(1:nq,indices_eig(indices_sort_ascent(1:nEF)))
    eigenvalues(1:nEF) = eigenvalues_temp(1:nEF)

    deallocate(work)
    deallocate(alpha_real)
    deallocate(alpha_imaginary)
    deallocate(vl)
    deallocate(vr)
    deallocate(beta)
    deallocate(eigenvalues_temp)
    deallocate(indices_eig_temp)
    deallocate(indices_eig)
    deallocate(indices_sort_ascent)
  end subroutine get_eigenvalues

  subroutine feast_message(info, boolean_abort)
    implicit none
    integer, intent(in) :: info
    logical, intent(inout) :: boolean_abort

    if (info .eq. 202) then
      print*, 'Error: Problem with size of the system n (n≤0)'
      boolean_abort = .TRUE.
    elseif (info .eq. 201) then
	    print*, 'Error: Problem with size of initial subspace m0 (m0≤0 or m0>n)'
      boolean_abort = .TRUE.
    elseif (info .eq. 200) then
      print*, 'Error: Problem with emin, emax (emin≥emax)'
      boolean_abort = .TRUE.
    elseif (info .eq. 100) then
      print*, 'Error: Problem with i-th value of the input Extended Eigensolver parameter (fpm(i)). Only the parameters in use are checked.'
      boolean_abort = .TRUE.
    elseif (info .eq. 4) then
      print*, 'Warning: Successful return of only the computed subspace after call with fpm(14)= 1'
      boolean_abort = .TRUE.
    elseif (info .eq. 3) then
      print*, 'Warning: Size of the subspace m0 is too small (m0<m)'
      boolean_abort = .TRUE.
   elseif (info .eq. 2) then
      print*, 'Warning: No Convergence (number of iteration loops >fpm(4))'
      boolean_abort = .TRUE.
    elseif (info .eq. 1) then
      print*, 'Warning: No eigenvalue found in the search interval. Change search intervall'
      boolean_abort = .TRUE.
    elseif (info .eq. 0) then
      print*, 'Succesful sparse calculation of eigenvalues'
    elseif (info .lt. 0) then
      print*, 'Error: Eigenvalue analysis failed!'
      boolean_abort = .TRUE.
     end if
  end subroutine feast_message

! count lines in an input file
  function count_lines(filename) result(nlines)
    character(LEN=*), intent(in) :: filename
    integer :: nlines, io

    open(99,file=filename)
    nlines = 0
    do
      read(99,*,iostat=io)
      if (io/=0) exit
      nlines = nlines + 1
    end do

    close(99)
  end function count_lines

! read file contents into array
 function readFileIntoArray(filename,nrows,ncolumns) result(array)
    character(LEN=*), intent(in) :: filename
    integer :: nrows, ncolumns, io_error, i_read
    real(kind = 8), allocatable, dimension(:,:) :: array

    allocate(array(nrows,ncolumns))

    open(unit = 110, file = filename, access='sequential', status = 'old', iostat = io_error)
      if (io_error == 0) then
        read(110, *)
        do i_read = 1, nrows
          read(unit=110, FMT='(10000000f30.15)') array(i_read,:)
        end do
      else
        print *, "ERROR: no file.dres file found!"
      end if
    close(110)

 end function readFileIntoArray

! write output in reduced form
 ! Subroutine for writing result files

  subroutine output_write_red(fileID,vector)
    real(kind = 8), dimension(:), allocatable, intent(in) :: vector
    integer :: fileID
3232 format(1000000f30.15)
    write(fileID,3232) vector(:)
  end subroutine output_write_red

  subroutine readprecalcfiles(q1, v1, lambda, oldsimfilename, boolean_IO_error)
    implicit none
    real(kind = 8), intent(inout), allocatable :: q1(:), v1(:), lambda(:)
    logical, intent(inout) :: boolean_IO_error
    character(:), allocatable :: oldsimfilename
    integer :: io_error

    if (allocated(oldsimfilename) .eqv. .FALSE.) then
      boolean_IO_error = .TRUE.
      return
    end if

    open(unit = 5000, file = oldsimfilename // '_qs.dres', status = 'old', iostat = io_error)
    if (io_error .eq. 0) then
      read(5000, *, iostat = io_error)
      do
        if (io_error .ne. 0) exit
        read(5000, *, iostat = io_error) q1(:)
      end do
      close(unit = 5000)
    else
      print*, '... error: no ', oldsimfilename // '_qs.dres', ' found!'
      boolean_IO_error = .TRUE.
    end if

    open(unit = 5000, file = oldsimfilename // '_vs.dres', status = 'old', iostat = io_error)
    if (io_error .eq. 0) then
      read(5000, *, iostat = io_error)
      do
        if (io_error .ne. 0) exit
        read(5000, *, iostat = io_error) v1(:)
      end do
      close(unit = 5000)
    else
      print*, '... warning: no ', oldsimfilename // '_vs.dres', ' found!'
    end if

    open(unit = 5000, file = oldsimfilename // '_lambdas.dres', status = 'old', iostat = io_error)
    if (io_error .eq. 0) then
      read(5000, *, iostat = io_error)
      do
        if (io_error .ne. 0) exit
        read(5000, *, iostat = io_error) lambda(:)
      end do
      close(unit = 5000)
    else
      print*, '... warning: no ', oldsimfilename // '_lambdas.dres', ' found!'
    end if

    return
  end subroutine readprecalcfiles

end module solver_functions
