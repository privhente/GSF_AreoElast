module solver_variables
  use class_sparse_matrix
  use sparse_function
  
  implicit none
  
  real(kind = 8), dimension(:), allocatable :: q0, q1, v0, v1, p1
  
  integer :: iteration, iterationlimit, istep, inos
  
  ! residuum and corresponding tolerance for the solution
  real(kind = 8) :: residuum !< the iteration residuum
  real(kind = 8), dimension(:), allocatable :: rhs_unbalance_force, &      !< residual vector for solution of iteration matrix
                                               lambda, lambda0, &                  !< lagrange multiplier vector
                                               deltaxvector                !< correction vector for the positional, velocity and lagrange multiplier solutions
  
  real(kind = 8), dimension(:), allocatable :: rhs_unbalance_force_sol
  real(kind = 8) :: rhs_unbalance_force_deviation, rhs_unbalance_force_amplitude
             
  type(sparse_matrix) :: sparse_smatrix             !< the sparse iteration matrix.
    
  integer, dimension(:), allocatable :: indicesq, &   !< positional (and director) indices into the iteration matrix
                                        indicesv, &   !< velocity (and director velocity) indices into the iteration matrix
                                        indicesg, &   !< lagrange multiplier equation indices into the iteration matrix
                                        indicesg_stat, &
                                        indicesq_mod  !< coordinates without shell strain dofs
  
! Arrays for extracting matrices
  real(kind = 8), dimension(:,:), allocatable   :: mass_matrix, stiffness_matrix, constraint_matrix, Hq_matrix, Hv_matrix
  real(kind = 8), dimension(:), allocatable :: sparse_csr_values_K, sparse_csr_values_M, sparse_csr_values_H, sparse_csr_values_H_t
  integer, dimension(:), allocatable :: sparse_csr_columns_K, sparse_csr_rowIndices_K, sparse_csr_compr_source_K
  integer, dimension(:), allocatable :: sparse_csr_columns_M, sparse_csr_rowIndices_M, sparse_csr_compr_source_M 
  integer, dimension(:), allocatable :: sparse_csr_columns_H, sparse_csr_rowIndices_H, sparse_csr_compr_source_H
  integer, dimension(:), allocatable :: sparse_csr_columns_H_t, sparse_csr_rowIndices_H_t, sparse_csr_compr_source_H_t

! general variables
  integer :: nq, nc, nn
  integer :: info
  integer :: convert_job(8)
    
  real(kind = 8), dimension(:), allocatable :: gravity_vector, &  !< gravity acceleration vector
                                              ml_t,            &  !< material load vector
                                              sl_t                !< spatial load vector
  
  integer :: gravityflag
  integer :: io_error_simulation
  
  real(kind = 8) :: deltat, totalt, time
  
  integer :: i, j, i6, i12, j6, j12, k, incoordinates
  
  integer :: tnnodes, tncoordinates, trconstraints
  
  integer :: row, col, index, nMatSize
    
  integer :: error
    
  integer, dimension(:), allocatable :: perm, perm1, perm2        !< permutation chosen by reordering
  
  real(kind = 8), allocatable :: temp_check(:)
  
end module solver_variables
