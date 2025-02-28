subroutine modal_solver(the_model_structure, settings, boolean_output_write)
! Importing classes
  use solver_functions
  use mkl_dss
  use class_model_structure, only: model_structure, step, append_int_vec
  use solver_variables
  use sparse_function
  use my_FileIO_functions
  
  implicit none

  type(model_structure), intent(inout) :: the_model_structure
  type(step), intent(inout)  :: settings
  type(single_sparse_matrix) :: sparse_nullspace

! Variables for solving eigenvalue problem
  integer :: nEF, nEF_update
  integer, allocatable :: indices_sort(:)
  real(kind = 8), allocatable, dimension(:)   :: eigenfrequencies, eigenvalues, eigenvalues_temp, dq
  real(kind = 8), allocatable, dimension(:,:) :: eigenvectors_red, eigenvectors, eigenvectors_red_temp
  real(kind = 8), dimension(:,:), allocatable :: Mred, Kred, N_matrix

! Arrays for sparse multiplication to get reduced matrices
  integer, allocatable :: sparse_csr_columns_N(:)

  integer, allocatable :: sparse_csr_columns_kred(:)
  integer, allocatable :: sparse_csr_columns_mred(:)
  integer, allocatable :: sparse_csr_columns_mred_temp(:)
  integer, allocatable :: sparse_csr_columns_kred_temp(:)
  integer, allocatable :: sparse_csr_columns_c(:)

  integer, allocatable :: sparse_csr_rowIndices_kred(:)
  integer, allocatable :: sparse_csr_rowIndices_mred(:)
  integer, allocatable :: sparse_csr_rowIndices_N(:)
  integer, allocatable :: sparse_csr_rowIndices_c(:)

  real(kind = 8), allocatable :: sparse_csr_values_kred(:)
  real(kind = 8), allocatable :: sparse_csr_values_mred(:)
  real(kind = 8), allocatable :: sparse_csr_values_kred_temp(:)
  real(kind = 8), allocatable :: sparse_csr_values_mred_temp(:)
  real(kind = 8), allocatable :: sparse_csr_values_N(:)
  real(kind = 8), allocatable :: sparse_csr_values_c(:)

! general variables
  integer :: num_nonzeros

! Indices arrays
  integer, dimension(:), allocatable :: indices_qg

! sparse solver variables
  integer :: loop, m, m0, int_eps_exp
  real(kind = 8) :: epsout, emin, emax
  real(kind = 8), allocatable :: res(:)
  integer :: fpm(128)
  integer, allocatable, dimension(:) :: indices_sort_ascent

  logical, intent(in) :: boolean_output_write
  
  convert_job(1) = 1
  convert_job(2) = 1
  convert_job(3) = 1
  convert_job(4) = 2

  tncoordinates  = the_model_structure%tncoordinates
  trconstraints  = the_model_structure%trconstraints

! Extracting only DOFs regarded to q, removing shell strain DOFs - 1. q rigid bodies, 2. q beams, 3. q shells, 4. velocity, 5. constraints
  nMatSize = the_model_structure%nnodes6*6 + the_model_structure%nnodes12*12 + trconstraints
  allocate(indices_qg(nMatSize))
  indices_qg = reshape((/indicesq_mod,indicesg/),(/nMatSize/))

! number of indices for coordinates and constrains
  nEF = settings%nEF
  nq  = the_model_structure%tncoordinates6 + the_model_structure%tncoordinates12
  nc  = trconstraints
  nn  = the_model_structure%tnnullspace

  print '(2X,A29)', "Starting modal calculation..."
  
! create sparse stiffness and mass matrices residual vector of last converged load step of dynamic calculation
  call sparse_smatrix%startAssembly()
  call the_model_structure%modelboundaries(the_model_structure%time_t)
  call the_model_structure%modelloads(the_model_structure%time_t, ml_t, sl_t)
  call the_model_structure%modelfk(the_model_structure%q_t, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, the_model_structure%q_t)

  call the_model_structure%modelgdg(the_model_structure%q_t, the_model_structure%q_t, sparse_smatrix, indicesg, indicesq, the_model_structure%lambda_t)

  ! Extracting stiffness and mass matrices from Iteration matrix
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, sparse_csr_values_K, sparse_csr_columns_K, sparse_csr_rowIndices_K, sparse_csr_compr_source_K, indicesq_mod, indicesq_mod)
  call extract_csr_Format_From_Sparsematrix(the_model_structure%sparse_mass, sparse_csr_values_M, sparse_csr_columns_M,sparse_csr_rowIndices_M, sparse_csr_compr_source_M, indicesq_mod,indicesq_mod)

! create sparse null space matrix
  print '(3X,A25,1X,I100)', "Calculating null space...", nn

  call sparse_initialize(sparse_nullspace,nq,nn)
  call the_model_structure%modelndgq(the_model_structure%q_t, the_model_structure%q_t,sparse_nullspace)
  call create_csr_format(sparse_nullspace)

  allocate(sparse_csr_values_N(size(sparse_nullspace%csr_values)))
  allocate(sparse_csr_columns_N(size(sparse_nullspace%csr_columns)))
  allocate(sparse_csr_rowIndices_N(size(sparse_nullspace%csr_rowIndices)))

  allocate(sparse_csr_values_c(nq*nn))
  allocate(sparse_csr_columns_c(nq*nn))
  allocate(sparse_csr_rowIndices_c(nn+1))

  allocate(sparse_csr_values_kred_temp(nn*nn))
  allocate(sparse_csr_columns_kred_temp(nn*nn))
  allocate(sparse_csr_rowIndices_kred(nn+1))

  allocate(sparse_csr_values_mred_temp(nn*nn))
  allocate(sparse_csr_columns_mred_temp(nn*nn))
  allocate(sparse_csr_rowIndices_mred(nn+1))

  sparse_csr_values_N     = sparse_nullspace%csr_values
  sparse_csr_columns_N    = sparse_nullspace%csr_columns
  sparse_csr_rowIndices_N = sparse_nullspace%csr_rowIndices

  print '(3X,A31)', "Calculating reduced matrices..."
  
! Performing sparse multiplication to get reduced stiffness matrix using mkl-function
  ! C = N_matrx' * S; C (nn x nq)
  sparse_csr_values_c = 0.0d0
  sparse_csr_columns_c = 0
  sparse_csr_rowIndices_c = 0
  call mkl_dcsrmultcsr('T', 0, 8, nq, nn, nq, sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, &
                                                      2.0d0*sparse_csr_values_K, sparse_csr_columns_K, sparse_csr_rowIndices_K, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, nq*nn, info)
! D = C*N_matrix; D (nn x nn)
  sparse_csr_values_kred_temp = 0.0d0
  sparse_csr_columns_kred_temp = 0
  sparse_csr_rowIndices_kred = 0
  call mkl_dcsrmultcsr('N', 0, 8, nn, nq, nn, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, &
                                                      sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, sparse_csr_values_kred_temp, sparse_csr_columns_kred_temp, sparse_csr_rowIndices_kred, nn*nn, info)

! Performing sparse multiplication to get reduced mass matrix using mkl-function
  ! C = N_matrx' * S; C (nn x nq)
  sparse_csr_values_c = 0.0d0
  sparse_csr_columns_c = 0
  sparse_csr_rowIndices_c = 0
  call mkl_dcsrmultcsr('T', 0, 8, nq, nn, nq, sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, &
                                                      sparse_csr_values_M, sparse_csr_columns_M, sparse_csr_rowIndices_M, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, nq*nn, info)
! D = C*N_matrix; D (nn x nn)
  sparse_csr_values_mred_temp = 0.0d0
  sparse_csr_columns_mred_temp = 0
  sparse_csr_rowIndices_mred = 0
  call mkl_dcsrmultcsr('N', 0, 8, nn, nq, nn, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, &
                                                      sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, sparse_csr_values_mred_temp, sparse_csr_columns_mred_temp, sparse_csr_rowIndices_mred, nn*nn, info)

  allocate(sparse_csr_values_kred(sparse_csr_rowIndices_kred(nn+1)-1))
  allocate(sparse_csr_values_mred(sparse_csr_rowIndices_mred(nn+1)-1))
  allocate(sparse_csr_columns_kred(sparse_csr_rowIndices_kred(nn+1)-1))
  allocate(sparse_csr_columns_mred(sparse_csr_rowIndices_mred(nn+1)-1))

  sparse_csr_values_kred = sparse_csr_values_kred_temp(1:sparse_csr_rowIndices_kred(nn+1)-1)
  sparse_csr_values_mred = sparse_csr_values_mred_temp(1:sparse_csr_rowIndices_mred(nn+1)-1)
  sparse_csr_columns_kred = sparse_csr_columns_kred_temp(1:sparse_csr_rowIndices_kred(nn+1)-1)
  sparse_csr_columns_mred = sparse_csr_columns_mred_temp(1:sparse_csr_rowIndices_mred(nn+1)-1)

! Calculating eigenfrequencies of augmented system
  num_nonzeros = max(size(sparse_csr_values_kred),size(sparse_csr_values_mred))
  nEF_update = nn

  print '(3X,A31)', "Calculating eigenfrequencies..."

  if (settings%sparse_flag == 1) then  ! Solving the eigenvalue problem using mkl - sparse - procedure
    print '(3X,A22)', "Using sparse solver..."
    m0           = nEF
    m            = 0
    emin         = (2.0d0*3.141592653589793d0)*(settings%emin**2)
    emax         = (2.0d0*3.141592653589793d0)*(settings%emax**2)
    int_eps_exp  = abs(int(log10(settings%tolerance)))

    allocate(eigenvalues_temp(m0))
    allocate(eigenvectors_red_temp(nn,m0))
    allocate(res(m0))
    res = 0.0d0

    call feastinit(fpm)
    !fpm(1) = 1
    fpm(2) = 16
    fpm(3) = int_eps_exp
    fpm(6) = 1

    call dfeast_scsrgv('F', nn, sparse_csr_values_kred, sparse_csr_rowIndices_kred, sparse_csr_columns_kred, &
                                sparse_csr_values_mred, sparse_csr_rowIndices_mred, sparse_csr_columns_mred, &
                                fpm, epsout, loop, emin, emax, m0, eigenvalues_temp, eigenvectors_red_temp, m, res, info)

    call feast_message(info, the_model_structure%boolean_abort)
    if (the_model_structure%boolean_abort .eqv. .FALSE.) then
      print '(3X,A36,1X,I100)', "Eigenvalue found in prescibed range:", m0
      nEF_update = m0
      call linSorting_real(indices_sort_ascent,eigenvalues_temp(1:nEF_update))
      allocate(eigenvalues(nEF_update))
      allocate(eigenvectors_red(nn,nEF_update))
      eigenvalues      = eigenvalues_temp(1:nEF_update)
      eigenvectors_red = 0.0d0
      do i=1,nEF_update
        eigenvectors_red(:,i) = eigenvectors_red_temp(:,indices_sort_ascent(i))
      end do
      deallocate(eigenvalues_temp)
      deallocate(eigenvectors_red_temp)
    end if
  else   ! Solving the eigenvalue problem using mkl - dense - procedure
    print '(3X,A21)', "Using dense solver..."
    allocate(Mred(nn,nn))
    allocate(Kred(nn,nn))

    Kred(:,:) = 0.0d0
    call mkl_ddnscsr(convert_job, nn, nn, Kred, nn, sparse_csr_values_kred, sparse_csr_columns_kred, sparse_csr_rowIndices_kred, info)

    Mred(:,:) = 0.0d0
    call mkl_ddnscsr(convert_job, nn, nn, Mred, nn, sparse_csr_values_mred, sparse_csr_columns_mred, sparse_csr_rowIndices_mred, info)

! Sort eigenvalues in ascending order
    call get_eigenvalues(eigenvalues, eigenvectors_red, Kred, Mred, nEF, nn, the_model_structure%boolean_abort)
    nEF_update = nEF
  end if

  if (the_model_structure%boolean_abort .neqv. .TRUE.) then
    allocate(N_matrix(nq,nn))
    call mkl_ddnscsr(convert_job, nq, nn, N_matrix, nq, sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, info)

    allocate(eigenvectors(nq,nEF_update))
    eigenvectors = 0.0d0
    if (nn == nq) then
      eigenvectors = eigenvectors_red
    else
      ! here inserting sparse multiplication
      eigenvectors = matmul(N_matrix,eigenvectors_red(:,1:nEF_update))
    end if

! Calculating eigenfrequencies: loop over eigenvalues to consider and signal negative eigenvalues in result file
    do i = 1, nEF_update
        if (eigenvalues(i) .gt. -1.0d0) then
            call append_int_vec(indices_sort,i)
        end if
    end do

    if (allocated(indices_sort) .eqv. .FALSE.) then
      the_model_structure%boolean_abort = .TRUE.
      print*,('Error: no eigenfrequencies found!')
      return
    end if
    
    nEF_update = size(indices_sort)
    allocate(eigenfrequencies(nEF_update))
    eigenfrequencies = 0.0d0
    do i = 1, nEF_update
      if (eigenvalues(indices_sort(i)) .lt. 0.0d0) then
        eigenfrequencies(i) = -sqrt(abs(eigenvalues(indices_sort(i))))/(2.0d0*3.141592653589793d0)
      else
        eigenfrequencies(i) = sqrt(abs(eigenvalues(indices_sort(i))))/(2.0d0*3.141592653589793d0)
      end if
    end do
    
    nEF = settings%nEF
    if (nEF > nEF_update) then
      nEF = nEF_update
    end if

! writing results of modal analysis
    print '(3X,A52)', "Writing results: eigenvector and eigenfrequencies..."
    print *, eigenfrequencies(1:nEF)
    allocate(dq(size(indicesq)))
    do i=1,nEF
      dq(:) = 0.0d0
      dq(indicesq_mod) = eigenvectors(:,indices_sort(i))

! writing in result file
      if (boolean_output_write .eqv. .TRUE.) then
        call the_model_structure%output_write(settings%step, eigenfrequencies(i), eigenfrequencies(i), i, settings%i_simutype, 0, dq, the_model_structure%qdot_t, the_model_structure%lambda_t, the_model_structure%stress)
      end if
    end do

! writing result for DeSiO test log file
    if (allocated(eigenvalues_temp)) then
      deallocate(eigenvalues_temp)
    end if
    allocate(eigenvalues_temp(1))
    eigenvalues_temp = eigenfrequencies(1)
    call writeRealVectorToFile(eigenvalues_temp,9999,'check.log')

    if (settings%output_flag == 1) then
      print *, "Output: Writing matrices..."

      call writeIntVectorToFile(indicesq_mod,   9999,'indicesq_mod.dres')
      call writeIntVectorToFile(indicesg,       9999,'indicesg.dres')

      call writeIntVectorToFile(sparse_csr_columns_M,   9999,'csr_columns_M.dres')
      call writeIntVectorToFile(sparse_csr_rowIndices_M,9999,'csr_rowIndices_M.dres')
      call writeRealVectorToFile(sparse_csr_values_M,   9999,'csr_values_M.dres')

      call writeIntVectorToFile(sparse_csr_columns_K,   9999,'csr_columns_K.dres')
      call writeIntVectorToFile(sparse_csr_rowIndices_K,9999,'csr_rowIndices_K.dres')
      call writeRealVectorToFile(sparse_csr_values_K,   9999,'csr_values_K.dres')

      call extract_csr_Format_From_Sparsematrix(sparse_smatrix, sparse_csr_values_H, sparse_csr_columns_H, sparse_csr_rowIndices_H, sparse_csr_compr_source_H, indicesg, indicesq_mod)
      call writeIntVectorToFile(sparse_csr_columns_H,   9999,'csr_columns_H.dres')
      call writeIntVectorToFile(sparse_csr_rowIndices_H,9999,'csr_rowIndices_H.dres')
      call writeRealVectorToFile(sparse_csr_values_H,   9999,'csr_values_H.dres')

      call writeIntVectorToFile(sparse_csr_columns_N,   9999,'csr_columns_N.dres')
      call writeIntVectorToFile(sparse_csr_rowIndices_N,9999,'csr_rowIndices_N.dres')
      call writeRealVectorToFile(sparse_csr_values_N,   9999,'csr_values_N.dres')
    end if

! Deallocating arrays
    deallocate(dq)
    if (allocated(Mred)) then
      deallocate(Mred)
      deallocate(Kred)
    end if
    
    deallocate(indices_sort)
    deallocate(N_matrix)
    deallocate(indices_qg)
    deallocate(eigenvalues)
    deallocate(eigenvectors)
    deallocate(eigenvectors_red)
    deallocate(sparse_csr_values_N)
    deallocate(sparse_csr_columns_N)
    deallocate(sparse_csr_rowIndices_N)
    deallocate(sparse_csr_values_c)
    deallocate(sparse_csr_columns_c)
    deallocate(sparse_csr_rowIndices_c)
    deallocate(sparse_csr_values_kred)
    deallocate(sparse_csr_columns_kred)
    deallocate(sparse_csr_rowIndices_kred)
    deallocate(sparse_csr_values_mred)
    deallocate(sparse_csr_columns_mred)
    deallocate(sparse_csr_rowIndices_mred)

    deallocate(sparse_csr_values_kred_temp)
    deallocate(sparse_csr_columns_kred_temp)
    deallocate(sparse_csr_values_mred_temp)
    deallocate(sparse_csr_columns_mred_temp)

  else

    if (allocated(eigenvalues_temp)) then
      deallocate(eigenvalues_temp)
    end if
    allocate(eigenvalues_temp(1))
    call writeRealVectorToFile(eigenvalues_temp*0.0d0,9999,'check.log')

  end if
  if (allocated(eigenfrequencies)) then
    deallocate(eigenfrequencies)
  end if

  return
  end subroutine modal_solver
