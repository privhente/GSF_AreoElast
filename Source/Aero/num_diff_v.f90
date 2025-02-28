subroutine num_diff_v(the_model_aero)
  use class_model_aero, only: model_aero
  use my_FileIO_functions, only: writeRealMatrixToFile
  implicit none
  
  type(model_aero), intent(inout) :: the_model_aero
  real(kind = 8), allocatable :: fae1(:), fae2(:), k_n(:,:), k_a(:,:), q0ae(:), v0ae(:), amatrix(:,:), oldcirculations(:), circulations(:), dfae(:,:)
  integer :: convert_job_csr_to_den(8), j, info
  real(kind = 8) :: epsil
  
  epsil = 1.0d-6
  
  convert_job_csr_to_den(1) = 1
  convert_job_csr_to_den(2) = 1
  convert_job_csr_to_den(3) = 1
  convert_job_csr_to_den(4) = 2
  
  allocate(q0ae(the_model_aero%tncoordinates))
  allocate(v0ae(the_model_aero%tncoordinates))
  allocate(amatrix(the_model_aero%tnrings, the_model_aero%tnrings))
  allocate(fae1(3*the_model_aero%tnrings))
  allocate(fae2(3*the_model_aero%tnrings))
  allocate(dfae(3*the_model_aero%tnrings,2*the_model_aero%tncoordinates))
  allocate(k_n(3*the_model_aero%tnrings,the_model_aero%tncoordinates))
  allocate(k_a(3*the_model_aero%tnrings,the_model_aero%tncoordinates))
  allocate(oldcirculations(the_model_aero%tnrings))
  allocate(circulations(the_model_aero%tnrings))
    
  q0ae            = the_model_aero%qs_t
  v0ae            = the_model_aero%vs_t
  circulations    = the_model_aero%circulations
  oldcirculations = the_model_aero%oldcirculations
  
!< calculating analytical tangent matrix
  the_model_aero%qs_t = q0ae
  the_model_aero%vs_t = v0ae
  the_model_aero%circulations = circulations
  the_model_aero%oldcirculations = oldcirculations
  
  call the_model_aero%actualize(.TRUE., .TRUE.)
  the_model_aero%airflow%vinfinity = the_model_aero%airflow%dir * the_model_aero%evaluateamplitude(the_model_aero%time)
  call the_model_aero%evaluaterhs(the_model_aero%cutoff)
  call the_model_aero%evaluateamatrix(the_model_aero%cutoff)
  the_model_aero%circulations = the_model_aero%rhs
  amatrix = the_model_aero%amatrix ! this has to be done here, when further using of amatrix
  call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
  the_model_aero%oldcirculations = oldcirculations
  call the_model_aero%actualize(.FALSE., .TRUE.)  ! actualize ring circulations
  call the_model_aero%modelfae()
  call the_model_aero%modeldfae(.TRUE.)
  
  call writeRealMatrixToFile(the_model_aero%dfae_ring,1000,'dfae_a.dres')
  
  !< calculating numerical tangent matrix  
  k_n = 0.0d0
  do j = 1,the_model_aero%tncoordinates
    the_model_aero%qs_t    = q0ae
    the_model_aero%vs_t    = v0ae
    the_model_aero%vs_t(j) = the_model_aero%vs_t(j) + epsil
    the_model_aero%circulations = circulations
    the_model_aero%oldcirculations = oldcirculations
    
    call the_model_aero%actualize(.TRUE., .TRUE.)
    the_model_aero%airflow%vinfinity = the_model_aero%airflow%dir * the_model_aero%evaluateamplitude(the_model_aero%time)
    call the_model_aero%evaluaterhs(the_model_aero%cutoff)
    call the_model_aero%evaluateamatrix(the_model_aero%cutoff)
    the_model_aero%circulations = the_model_aero%rhs
    amatrix = the_model_aero%amatrix ! this has to be done here, when further using of amatrix
    call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
    the_model_aero%oldcirculations = oldcirculations
    call the_model_aero%actualize(.FALSE., .TRUE.)
    call the_model_aero%modelfae()
    fae1 = the_model_aero%fae
    
    the_model_aero%qs_t    = q0ae
    the_model_aero%vs_t    = v0ae
    the_model_aero%vs_t(j) = the_model_aero%vs_t(j) - epsil
    the_model_aero%circulations = circulations
    the_model_aero%oldcirculations = oldcirculations
    
    call the_model_aero%actualize(.TRUE., .TRUE.)
    the_model_aero%airflow%vinfinity = the_model_aero%airflow%dir * the_model_aero%evaluateamplitude(the_model_aero%time)
    call the_model_aero%evaluaterhs(the_model_aero%cutoff)
    call the_model_aero%evaluateamatrix(the_model_aero%cutoff)
    the_model_aero%circulations = the_model_aero%rhs
    amatrix = the_model_aero%amatrix ! this has to be done here, when further using of amatrix
    call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
    the_model_aero%oldcirculations = oldcirculations
    call the_model_aero%actualize(.FALSE., .TRUE.)
    call the_model_aero%modelfae()
    fae2 = the_model_aero%fae
   
    k_n(1:3*the_model_aero%tnrings,j) = ( fae1 - fae2 ) / (2.0d0*epsil)   
  end do
  
  call writeRealMatrixToFile(k_n,1000,'dfaev_n.dres')

  print*, maxval(abs(the_model_aero%dfae_ring(1:3*the_model_aero%tnrings,the_model_aero%tncoordinates+1:2*the_model_aero%tncoordinates)-k_n))
  return
end subroutine num_diff_v