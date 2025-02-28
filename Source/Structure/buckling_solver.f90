subroutine buckling_solver(the_model_structure, settings, boolean_output_write)
! Declaring variables
  use solver_functions

  use class_model_structure, only: model_structure, step, append_int_vec
  use solver_variables
  use sparse_function
  use my_FileIO_functions
  
  implicit none

  type(model_structure), intent(inout) :: the_model_structure
  type(step), intent(inout)  :: settings
  type(single_sparse_matrix) :: sparse_nullspace           !< the sparse iteration matrix.

  logical, intent(in) :: boolean_output_write

! Variables for solving eigenvalue problem
  integer :: nEF, nEF_update
  integer, allocatable :: indices_sort(:)
  real(kind = 8), allocatable, dimension(:)   :: eigenvalues, eigenvalues_temp, dq
  real(kind = 8), allocatable, dimension(:,:) :: eigenvectors_red, eigenvectors, eigenvectors_red_temp
  real(kind = 8), dimension(:,:), allocatable :: Kmat_red, Kgeo_red, N_matrix

! Arrays for sparse multiplication to get reduced matrices
  integer, allocatable :: sparse_csr_columns_N(:)
  integer, allocatable :: sparse_csr_columns_c(:)
  integer, allocatable :: sparse_csr_columns_Kmat(:)
  integer, allocatable :: sparse_csr_columns_Kmat_q(:)
  integer, allocatable :: sparse_csr_columns_Kgeo(:)
  integer, allocatable :: sparse_csr_columns_kmat_red(:)
  integer, allocatable :: sparse_csr_columns_kgeo_red(:)

  integer, allocatable :: sparse_csr_columns_kmat_red_temp(:)
  integer, allocatable :: sparse_csr_columns_kgeo_red_temp(:)

  integer, allocatable :: sparse_csr_rowIndices_N(:)
  integer, allocatable :: sparse_csr_rowIndices_c(:)
  integer, allocatable :: sparse_csr_rowIndices_Kmat(:)
  integer, allocatable :: sparse_csr_rowIndices_Kmat_q(:)
  integer, allocatable :: sparse_csr_rowIndices_Kgeo(:)
  integer, allocatable :: sparse_csr_rowIndices_kmat_red(:)
  integer, allocatable :: sparse_csr_rowIndices_kgeo_red(:)

  integer, allocatable :: sparse_csr_compr_source_Kmat(:)
  integer, allocatable :: sparse_csr_compr_source_Kmat_q(:)
  integer, allocatable :: sparse_csr_compr_source_Kgeo(:)

  real(kind = 8), allocatable :: sparse_csr_values_N(:)
  real(kind = 8), allocatable :: sparse_csr_values_c(:)
  real(kind = 8), allocatable :: sparse_csr_values_kmat(:)
  real(kind = 8), allocatable :: sparse_csr_values_kmat_q(:)
  real(kind = 8), allocatable :: sparse_csr_values_kgeo(:)
  real(kind = 8), allocatable :: sparse_csr_values_kmat_red(:)
  real(kind = 8), allocatable :: sparse_csr_values_kgeo_red(:)
  real(kind = 8), allocatable :: sparse_csr_values_kmat_red_temp(:)
  real(kind = 8), allocatable :: sparse_csr_values_kgeo_red_temp(:)

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

  print '(2X,A32)', "Starting buckling calculation..."

! Material stiffness matrix
  call the_model_structure%localSimuSetting(1.0d0,2,.FALSE.,.TRUE.)         ! the_model_structure%localSimuSetting(delta, simutype, kgeo_on, kmat_on)
  call sparse_smatrix%startAssembly()
  call the_model_structure%modelloads(0.0d0, ml_t, sl_t)
  call the_model_structure%modelfk(the_model_structure%q_0, ml_t, sl_t, sparse_smatrix, indicesq, indicesv)
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, sparse_csr_values_kmat, sparse_csr_columns_kmat, sparse_csr_rowIndices_kmat, sparse_csr_compr_source_kmat, indicesq_mod, indicesq_mod)

! Kgeo aus Kmat(q)
  call the_model_structure%localSimuSetting(1.0d0,2,.FALSE.,.TRUE.)         ! the_model_structure%localSimuSetting(delta, simutype, kgeo_on, kmat_on)
  call sparse_smatrix%startAssembly()
  call the_model_structure%modelloads(0.0d0, ml_t, sl_t)
  call the_model_structure%modelfk(the_model_structure%q_t, ml_t, sl_t, sparse_smatrix, indicesq, indicesv)
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, sparse_csr_values_kmat_q, sparse_csr_columns_kmat_q, sparse_csr_rowIndices_kmat_q, sparse_csr_compr_source_kmat_q, indicesq_mod, indicesq_mod)
!  sparse_csr_values_kmat_q = sparse_csr_values_kmat - sparse_csr_values_kmat_q

! Geometrical stiffness matrix and Hessian matrix of constraints
  call the_model_structure%localSimuSetting(1.0d0,2,.TRUE.,.FALSE.)         ! the_model_structure%localSimuSetting(delta, simutype, kgeo_on, kmat_on)
  call sparse_smatrix%startAssembly()
  call the_model_structure%modelboundaries(the_model_structure%time_t)
  call the_model_structure%modelloads(the_model_structure%time_t, ml_t, sl_t)
  call the_model_structure%modelfk(the_model_structure%q_t, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, the_model_structure%q_t)
  call the_model_structure%modelgdg(the_model_structure%q_t, the_model_structure%q_t, sparse_smatrix, indicesg, indicesq, the_model_structure%lambda_t)
  call extract_csr_Format_From_Sparsematrix(sparse_smatrix, sparse_csr_values_kgeo, sparse_csr_columns_kgeo, sparse_csr_rowIndices_kgeo, sparse_csr_compr_source_kgeo, indicesq_mod, indicesq_mod)
  sparse_csr_values_kgeo = sparse_csr_values_kgeo !+ sparse_csr_values_kmat_q

  print '(3X,A25,1X,I100)', "Calculating null space...", nn
  
! create sparse null space matrix
  call sparse_initialize(sparse_nullspace,nq,nn)
  call the_model_structure%modelndgq(the_model_structure%q_t, the_model_structure%q_t,sparse_nullspace)
  call create_csr_format(sparse_nullspace)

  allocate(sparse_csr_values_N(size(sparse_nullspace%csr_values)))
  allocate(sparse_csr_columns_N(size(sparse_nullspace%csr_columns)))
  allocate(sparse_csr_rowIndices_N(size(sparse_nullspace%csr_rowIndices)))

! null-space-projection of the stiffness and mass matrices
  allocate(sparse_csr_values_c(nq*nn))
  allocate(sparse_csr_columns_c(nq*nn))
  allocate(sparse_csr_rowIndices_c(nn+1))

  allocate(sparse_csr_values_kmat_red_temp(nn*nn))
  allocate(sparse_csr_columns_kmat_red_temp(nn*nn))
  allocate(sparse_csr_rowIndices_kmat_red(nn+1))

  allocate(sparse_csr_values_kgeo_red_temp(nn*nn))
  allocate(sparse_csr_columns_kgeo_red_temp(nn*nn))
  allocate(sparse_csr_rowIndices_kgeo_red(nn+1))

  sparse_csr_values_N     = sparse_nullspace%csr_values
  sparse_csr_columns_N    = sparse_nullspace%csr_columns
  sparse_csr_rowIndices_N = sparse_nullspace%csr_rowIndices

  print '(3X,A31)', "Calculating reduced matrices..."

! Performing sparse multiplication to get reduced stiffness matrix using mkl-function
! C = N_matrx' * S; C (nn x nq)
  sparse_csr_values_c     = 0.0d0
  sparse_csr_columns_c    = 0
  sparse_csr_rowIndices_c = 0
  call mkl_dcsrmultcsr('T', 0, 8, nq, nn, nq, sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, &
                                                      sparse_csr_values_kmat, sparse_csr_columns_Kmat, sparse_csr_rowIndices_Kmat, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, nq*nn, info)
! D = C*N_matrix; D (nn x nn)
  sparse_csr_values_kmat_red_temp  = 0.0d0
  sparse_csr_columns_kmat_red_temp = 0
  sparse_csr_rowIndices_kmat_red   = 0
  call mkl_dcsrmultcsr('N', 0, 8, nn, nq, nn, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, &
                                                      sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, sparse_csr_values_kmat_red_temp, sparse_csr_columns_kmat_red_temp, sparse_csr_rowIndices_kmat_red, nn*nn, info)
! Performing sparse multiplication to get reduced mass matrix using mkl-function
  ! C = N_matrx' * S; C (nn x nq)
  sparse_csr_values_c     = 0.0d0
  sparse_csr_columns_c    = 0
  sparse_csr_rowIndices_c = 0
  call mkl_dcsrmultcsr('T', 0, 8, nq, nn, nq, sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, &
                                                      sparse_csr_values_Kgeo, sparse_csr_columns_Kgeo, sparse_csr_rowIndices_Kgeo, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, nq*nn, info)
! D = C*N_matrix; D (nn x nn)
  sparse_csr_values_kgeo_red_temp  = 0.0d0
  sparse_csr_columns_kgeo_red_temp = 0
  sparse_csr_rowIndices_kgeo_red   = 0
  call mkl_dcsrmultcsr('N', 0, 8, nn, nq, nn, sparse_csr_values_c, sparse_csr_columns_c, sparse_csr_rowIndices_c, &
                                                      sparse_csr_values_N, sparse_csr_columns_N, sparse_csr_rowIndices_N, sparse_csr_values_kgeo_red_temp, sparse_csr_columns_kgeo_red_temp, sparse_csr_rowIndices_kgeo_red, nn*nn, info)

  allocate(sparse_csr_values_kmat_red(sparse_csr_rowIndices_kmat_red(nn+1)-1))
  allocate(sparse_csr_values_kgeo_red(sparse_csr_rowIndices_kgeo_red(nn+1)-1))
  allocate(sparse_csr_columns_kmat_red(sparse_csr_rowIndices_kmat_red(nn+1)-1))
  allocate(sparse_csr_columns_kgeo_red(sparse_csr_rowIndices_kgeo_red(nn+1)-1))

  sparse_csr_values_kmat_red  = sparse_csr_values_kmat_red_temp(1:sparse_csr_rowIndices_kmat_red(nn+1)-1)
  sparse_csr_values_kgeo_red  = sparse_csr_values_kgeo_red_temp(1:sparse_csr_rowIndices_kgeo_red(nn+1)-1)
  sparse_csr_columns_kmat_red = sparse_csr_columns_kmat_red_temp(1:sparse_csr_rowIndices_kmat_red(nn+1)-1)
  sparse_csr_columns_kgeo_red = sparse_csr_columns_kgeo_red_temp(1:sparse_csr_rowIndices_kgeo_red(nn+1)-1)

! Calculating eigenfrequencies of augmented system, output in ascending order
  print '(3X,A27)', "Calculating load factors..."
  
  num_nonzeros = max(size(sparse_csr_values_kgeo_red),size(sparse_csr_values_kmat_red))
  nEF_update = nn
  if (settings%sparse_flag == 1) then  ! Solving the eigenvalue problem using mkl - sparse
    print '(3X,A22)', "Using sparse solver..."
    m0 = nEF
    m  = 0
    emin = settings%emin
    emax = settings%emax
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

    call dfeast_scsrgv('F', nn, sparse_csr_values_kmat_red, sparse_csr_rowIndices_kmat_red, sparse_csr_columns_kmat_red, &
                                -sparse_csr_values_kgeo_red, sparse_csr_rowIndices_kgeo_red, sparse_csr_columns_kgeo_red, &
                                                          fpm, epsout, loop, emin, emax, m0, eigenvalues_temp, eigenvectors_red_temp, m, res, info)
    call feast_message(info, the_model_structure%boolean_abort)
    if (the_model_structure%boolean_abort .eqv. .FALSE.) then
       nEF_update = m0
      call linSorting_real(indices_sort_ascent,eigenvalues_temp(1:nEF_update))
      allocate(eigenvalues(nEF_update))
      allocate(eigenvectors_red(nn,nEF_update))
      eigenvalues = eigenvalues_temp(1:nEF_update)
      do i=1,nEF_update
        eigenvectors_red(:,i) = eigenvectors_red_temp(:,indices_sort_ascent(i))
      end do
      deallocate(eigenvalues_temp)
      deallocate(eigenvectors_red_temp)
    end if

  else   ! Solving the eigenvalue problem using dense - procedure
    print '(3X,A21)', "Using dense solver..."
    allocate(Kmat_red(nn,nn))
    allocate(Kgeo_red(nn,nn))

    Kmat_red(:,:) = 0.0d0
    call mkl_ddnscsr(convert_job, nn, nn, Kmat_red, nn, sparse_csr_values_kmat_red, sparse_csr_columns_kmat_red, sparse_csr_rowIndices_kmat_red, info)

    Kgeo_red(:,:) = 0.0d0
    call mkl_ddnscsr(convert_job, nn, nn, Kgeo_red, nn, sparse_csr_values_kgeo_red, sparse_csr_columns_kgeo_red, sparse_csr_rowIndices_kgeo_red, info)

    call get_eigenvalues(eigenvalues,eigenvectors_red,Kmat_red,-Kgeo_red,nEF,nn,the_model_structure%boolean_abort)
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

  ! using only eigenvalues greater equal zero
      do i = 1, nEF_update
        if (eigenvalues(i) >= 0.0d0) then
          call append_int_vec(indices_sort, i)
        end if
      end do

      if (allocated(indices_sort) .eqv. .FALSE.) then
        the_model_structure%boolean_abort = .TRUE.
        print*,('Error: no eigenvalues found!')
        return
      end if
      
      nEF_update = size(indices_sort)
      nEF = settings%nEF
      if (nEF > nEF_update) then
        nEF = nEF_update
      end if

  ! writing results of linear buckling analysis
      print '(3X,A59)', "Writing results: buckling modes and critical load factor..."
      print *, eigenvalues(indices_sort(1:nEF))

      allocate(dq(size(indicesq)))
      do i = 1,nEF
        dq(:) = 0.0d0
        dq(indicesq_mod) = eigenvectors(:,indices_sort(i))
  ! writing in result file
        if (boolean_output_write .eqv. .TRUE.) then
          call the_model_structure%output_write(settings%step, eigenvalues(indices_sort(i)), eigenvalues(indices_sort(i)), i, settings%i_simutype, 0, dq, the_model_structure%qdot_t, the_model_structure%lambda_t, the_model_structure%stress)
        end if 
      end do

! writing result for DeSiO test log file
      if (allocated(eigenvalues_temp)) then
        deallocate(eigenvalues_temp)
      end if
      allocate(eigenvalues_temp(1))
      eigenvalues_temp = eigenvalues(indices_sort(1))
      call writeRealVectorToFile(eigenvalues_temp,9999,'check.log')

      if (settings%output_flag == 1) then
          print *, "Output: Writing matrices..."

          call writeIntVectorToFile(indicesq_mod,   9999,'indicesq_mod.dres')
          call writeIntVectorToFile(indicesg,       9999,'indicesg.dres')

          call writeIntVectorToFile(sparse_csr_columns_Kgeo,   9999,'csr_columns_Kgeo.dres')
          call writeIntVectorToFile(sparse_csr_rowIndices_Kgeo,9999,'csr_rowIndices_Kgeo.dres')
          call writeRealVectorToFile(sparse_csr_values_Kgeo,   9999,'csr_values_Kgeo.dres')

          call writeIntVectorToFile(sparse_csr_columns_Kmat,   9999,'csr_columns_Kmat.dres')
          call writeIntVectorToFile(sparse_csr_rowIndices_Kmat,9999,'csr_rowIndices_Kmat.dres')
          call writeRealVectorToFile(sparse_csr_values_kmat,   9999,'csr_values_Kmat.dres')

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
      if (allocated(Kgeo_red)) then
        deallocate(Kgeo_red)
        deallocate(Kmat_red)
      end if
      deallocate(N_matrix)
      deallocate(indices_qg)
      if (allocated(indices_sort)) then
        deallocate(indices_sort)
      end if
      deallocate(eigenvectors)
      deallocate(eigenvectors_red)
      deallocate(sparse_csr_values_N)
      deallocate(sparse_csr_columns_N)
      deallocate(sparse_csr_rowIndices_N)
      deallocate(sparse_csr_values_c)
      deallocate(sparse_csr_columns_c)
      deallocate(sparse_csr_rowIndices_c)
      deallocate(sparse_csr_values_kmat_red)
      deallocate(sparse_csr_columns_kmat_red)
      deallocate(sparse_csr_rowIndices_kmat_red)
      deallocate(sparse_csr_values_kgeo_red)
      deallocate(sparse_csr_columns_kgeo_red)
      deallocate(sparse_csr_rowIndices_kgeo_red)

      deallocate(sparse_csr_values_kmat_red_temp)
      deallocate(sparse_csr_columns_kmat_red_temp)
      deallocate(sparse_csr_values_kgeo_red_temp)
      deallocate(sparse_csr_columns_kgeo_red_temp)

  else

      if (allocated(eigenvalues_temp)) then
        deallocate(eigenvalues_temp)
      end if
      allocate(eigenvalues_temp(1))
      call writeRealVectorToFile(eigenvalues_temp*0.0d0,9999,'check.log')

  end if
  if (allocated(eigenvalues)) then
    deallocate(eigenvalues)
  end if

  return
  end subroutine buckling_solver
