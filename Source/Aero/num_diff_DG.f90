subroutine num_diff_DG(the_model_aero)
  use class_model_aero,  only: model_aero
  use my_FileIO_functions, only: writeRealMatrixToFile, writeRealVectorToFile
  implicit none
  
  type(model_aero), intent(inout) :: the_model_aero
  real(kind = 8), allocatable :: G1(:), G2(:), G_a(:,:), G_n(:,:)
  real(kind = 8), allocatable :: q0ae(:), v0ae(:), amatrix(:,:), oldcirculations(:), circulations(:)
  integer :: convert_job(8), j, info
  real(kind = 8) :: epsil
  
  real(kind = 8) :: dbl_a
  
  epsil = 1.0d-6
  convert_job(1) = 1
  convert_job(2) = 1
  convert_job(3) = 1
  convert_job(4) = 2

  allocate(q0ae(the_model_aero%tncoordinates))
  allocate(v0ae(the_model_aero%tncoordinates))
  allocate(G1(the_model_aero%tnrings))
  allocate(G2(the_model_aero%tnrings))
  allocate(G_n(the_model_aero%tnrings, the_model_aero%tncoordinates))
  allocate(G_a(the_model_aero%tnrings, the_model_aero%tncoordinates))
  allocate(amatrix(the_model_aero%tnrings, the_model_aero%tnrings))
  allocate(oldcirculations(the_model_aero%tnrings))
  allocate(circulations(the_model_aero%tnrings))
    
  q0ae            = the_model_aero%qs_t
  v0ae            = the_model_aero%vs_t
  circulations    = the_model_aero%circulations
  oldcirculations = the_model_aero%oldcirculations

!< calculating tangent matrix analytically
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
  call the_model_aero%actualize(.FALSE., .TRUE.)  ! actualize ring circulations
  call the_model_aero%modelfae()
  call the_model_aero%modeldfae(.TRUE.)
  call writeRealMatrixToFile(the_model_aero%dG(1:the_model_aero%tnrings,1:2*the_model_aero%tncoordinates),1000,'dG_a.dres')

!< calculating tangent matrix with respect to velocities numerically
  G_n = 0.0d0
  do j = 1,the_model_aero%tncoordinates
    the_model_aero%qs_t    = q0ae
    the_model_aero%qs_t(j) = the_model_aero%qs_t(j) + epsil
    the_model_aero%vs_t    = v0ae
    the_model_aero%circulations = circulations
    the_model_aero%oldcirculations = oldcirculations
    
    call the_model_aero%actualize(.TRUE., .TRUE.)
    the_model_aero%airflow%vinfinity = the_model_aero%airflow%dir * the_model_aero%evaluateamplitude(the_model_aero%time)
    call the_model_aero%evaluaterhs(the_model_aero%cutoff)
    call the_model_aero%evaluateamatrix(the_model_aero%cutoff)
    the_model_aero%circulations = the_model_aero%rhs
    amatrix = the_model_aero%amatrix
    call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
    G1 = the_model_aero%circulations
    call the_model_aero%actualize(.FALSE., .TRUE.)
    call the_model_aero%modelfae()
    
    the_model_aero%qs_t    = q0ae
    the_model_aero%qs_t(j) = the_model_aero%qs_t(j) - epsil
    the_model_aero%vs_t    = v0ae
    the_model_aero%circulations = circulations
    the_model_aero%oldcirculations = oldcirculations
    
    call the_model_aero%actualize(.TRUE., .TRUE.)
    the_model_aero%airflow%vinfinity = the_model_aero%airflow%dir * the_model_aero%evaluateamplitude(the_model_aero%time)
    call the_model_aero%evaluaterhs(the_model_aero%cutoff)
    call the_model_aero%evaluateamatrix(the_model_aero%cutoff)
    the_model_aero%circulations = the_model_aero%rhs
    amatrix = the_model_aero%amatrix
    call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
    G2 = the_model_aero%circulations
    call the_model_aero%actualize(.FALSE., .TRUE.)
    call the_model_aero%modelfae()
    
    G_n(1:the_model_aero%tnrings,j) =  (G1 - G2) / (2.0d0*epsil)
  end do
  call writeRealMatrixToFile(G_n,1000,'dGx_n.dres')
  print*, 'max diff dG_n = ', maxval(abs(the_model_aero%dG(1:the_model_aero%tnrings,1:the_model_aero%tncoordinates)-G_n))

!< calculating tangent matrix with respect to velocities numerically
  G_n = 0.0d0
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
    amatrix = the_model_aero%amatrix
    call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
    G1 = the_model_aero%circulations
    
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
    amatrix = the_model_aero%amatrix
    call dgesv(the_model_aero%tnrings, 1, amatrix, the_model_aero%tnrings, the_model_aero%ipiv, the_model_aero%circulations, the_model_aero%tnrings, info)
    G2 = the_model_aero%circulations
    G_n(1:the_model_aero%tnrings,j) =  (G1 - G2) / (2.0d0*epsil)
  end do
  call writeRealMatrixToFile(G_n,1000,'dGv_n.dres')
  print*, 'max diff dG_v = ', maxval(abs(the_model_aero%dG(1:the_model_aero%tnrings,the_model_aero%tncoordinates+1:2*the_model_aero%tncoordinates)-G_n))
  
  return
end subroutine num_diff_DG