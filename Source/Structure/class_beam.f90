module class_beam

  use class_node_12
  use class_beam_element
  use class_sparse_matrix

  implicit none
  
  type :: beam
     
     integer :: nnodes
     integer :: nelements
     integer :: ncoordinates
     type(node_12), allocatable :: nodes(:)
     type(beam_element), allocatable :: elements(:)
     real(kind = 8) :: penergy, deltat
     integer, allocatable :: indicesq12(:), indicesst(:)
     integer :: simutype
     real(kind = 8), dimension(:), allocatable :: f, fqint, fqdyn, fv, fqext, stress
     logical :: flag_kgeo_on, flag_kmat_on
   contains

     procedure :: initialization => beaminitialization
     procedure :: actualization => beamactualization
     procedure, private :: beamm
     procedure, private :: beamefk
     
  end type beam

contains

  subroutine beaminitialization(this, sparse_mass, offset)

    class(beam), intent(inout) :: this
    class(sparse_matrix), intent(inout) :: sparse_mass

    integer, intent(in) :: offset
    integer :: i, i12, j, j12
    
    this%ncoordinates = 12*this%nnodes

    allocate(this%fqint  (this%ncoordinates))
    allocate(this%fqdyn  (this%ncoordinates))
    allocate(this%fv     (this%ncoordinates))
    allocate(this%fqext  (this%ncoordinates))
    allocate(this%stress (6*this%nelements))
    
    do i = 1, this%nnodes

       i12 = 12*(i-1)
     
       this%nodes(i)%coordinates(:) = [i12+ 1, i12+ 2, i12+ 3, i12+ 4, i12+ 5, i12+ 6, i12+ 7, i12+ 8, i12+ 9, i12+10, i12+11, i12+12]
       
    end do
    
    do i = 1, this%nelements
       
       do j = 1, 2 
          
          j12 = 12*(j-1)
          
          this%elements(i)%q_0(j12+ 1:j12+12) = this%nodes(this%elements(i)%connectivity(j))%q_0
          
       end do

       this%elements(i)%deltat = this%deltat
       
    end do
        
    call this%beamm(sparse_mass, offset)

    return
    
  end subroutine beaminitialization

  subroutine beamactualization(this, q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2, v_1, v_2)

    class(beam), intent(inout) :: this
    real(kind = 8), intent(in) :: q_1(:)
    real(kind = 8), intent(in) :: ml_a(:)
    real(kind = 8), intent(in) :: sl_a(:)
    
    integer, intent(in) :: qOffset, vOffset
    type(sparse_matrix),  intent(inout)  :: sparse_smatrix
    
    real(kind = 8), optional, intent(in) :: q_2(:)
    real(kind = 8), optional, intent(in) :: v_1(:)
    real(kind = 8), optional, intent(in) :: v_2(:)

    if (present(q_2)) then

     if (present(v_1) .and. present(v_2)) then
        
        call this%beamefk(q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2, v_1, v_2)
        
     else
       
        call this%beamefk(q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2)

     end if

    else

       call this%beamefk(q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset)
       
    end if
    
  end subroutine beamactualization

  subroutine beamm(this, sparse_mass, offset)

    implicit none
    
    class(beam), intent(inout) :: this
    class(sparse_matrix), intent(inout) :: sparse_mass
    integer, intent(in) :: offset
    
    integer :: i, indices24(24)
    
    do i = 1, size(this%elements)
       
       indices24 = [this%nodes(this%elements(i)%connectivity(1))%coordinates(:), &
                    this%nodes(this%elements(i)%connectivity(2))%coordinates(:)]
              
       call sparse_mass%assemblyIndexed(indices24 + offset, indices24 + offset, this%elements(i)%massmatrix())
       
    end do
    
  end subroutine beamm
  
  subroutine beamefk(this, q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2, v_1, v_2)

    implicit none

    class(beam), intent(inout) :: this
    real(kind = 8), intent(in) :: q_1(:)
    real(kind = 8), intent(in) :: ml_a(:)
    real(kind = 8), intent(in) :: sl_a(:)
    
    integer, intent(in) :: qOffset, vOffset
    type(sparse_matrix),  intent(inout)  :: sparse_smatrix
    
    real(kind = 8), optional, intent(in) :: q_2(:)
    real(kind = 8), optional, intent(in) :: v_1(:)
    real(kind = 8), optional, intent(in) :: v_2(:)

    integer :: i, i6, indices6(6), indices12(12), indices24(24)
    real(kind = 8), dimension(24) :: temp_q1, temp_q2, temp_v1, temp_v2
    real(kind = 8), dimension(6) :: temp_ml_a, temp_sl_a
    real(kind = 8), dimension(12) :: temp_q1_12

    this%penergy   = 0.0d0
    this%fqint(:)  = 0.0d0
    this%fqdyn(:)  = 0.0d0
    this%fv(:)     = 0.0d0    
    this%stress(:) = 0.0d0
    
! evaluation of internal terms at the elements, elastic contribution.
    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(indices24)
    
    do i = 1, this%nelements
       indices24 = [this%nodes(this%elements(i)%connectivity(1))%coordinates, this%nodes(this%elements(i)%connectivity(2))%coordinates]
       this%elements(i)%deltat = this%deltat
       
       if (present(q_2)) then
          if (present(v_1) .and. present(v_2)) then
            temp_q1=q_1(indices24)
            temp_q2=q_2(indices24)
            temp_v1=v_1(indices24)
            temp_v2=v_2(indices24)
            !call this%elements(i)%internalterms(q_1(indices24), q_2(indices24), v_1(indices24), v_2(indices24))
            call this%elements(i)%internalterms(temp_q1, temp_q2, temp_v1, temp_v2)
          else
            temp_q1=q_1(indices24)
            temp_q2=q_2(indices24)
            !call this%elements(i)%internalterms(q_1(indices24), q_2(indices24))
            call this%elements(i)%internalterms(temp_q1, temp_q2)
          end if
       else
          call this%elements(i)%internalterms(q_1(indices24))
       end if
    end do

!$OMP  END PARALLEL DO
    
    do i = 1, this%nelements
      indices24 = [this%nodes(this%elements(i)%connectivity(1))%coordinates, this%nodes(this%elements(i)%connectivity(2))%coordinates]
! this is formulated in terms of addition, so there will be data corruption if run in parallel
      this%penergy           = this%penergy           + this%elements(i)%penergy
      this%fqint (indices24) = this%fqint (indices24) + this%elements(i)%fqint
      this%fqdyn (indices24) = this%fqdyn (indices24) + this%elements(i)%fqdyn
      this%fv (indices24)    = this%fv (indices24)    + this%elements(i)%fv

      this%stress(6*(i-1)+1:6*(i-1)+6) = this%elements(i)%stress_resultants
      
      call sparse_smatrix%assemblyIndexed(indices24 + qOffset, indices24 + qOffset, this%elements(i)%kqq)
      call sparse_smatrix%assemblyIndexed(indices24 + vOffset, indices24 + vOffset, this%elements(i)%kvv)
      call sparse_smatrix%assemblyIndexed(indices24 + qOffset, indices24 + vOffset, this%elements(i)%kqv)
      call sparse_smatrix%assemblyIndexed(indices24 + vOffset, indices24 + qOffset, this%elements(i)%kvq)
      
    end do
    
! evaluation of external terms at the nodes, material and spatial loads.

    this%fqext (:)    = 0.0d0
    
    do i = 1, this%nnodes

       i6 = 6*(i-1)
       
       indices6 = [i6+1, i6+2, i6+3, i6+4, i6+5, i6+6]
       
       indices12 = this%nodes(i)%coordinates
       
       temp_ml_a=ml_a(indices6)
       temp_sl_a=sl_a(indices6)

       if (present(q_2)) then
          !call this%nodes(i)%externalterms(0.5d0*(q_1(indices12)+q_2(indices12)), ml_a(indices6), sl_a(indices6))
          temp_q1_12=0.5d0*(q_1(indices12)+q_2(indices12))
          call this%nodes(i)%externalterms(temp_q1_12, temp_ml_a, temp_sl_a)
       else
          !call this%nodes(i)%externalterms(q_1(indices12), ml_a(indices6), sl_a(indices6))
          temp_q1_12=q_1(indices12)
          call this%nodes(i)%externalterms(temp_q1_12, temp_ml_a, temp_sl_a)
       end if

       this%fqext (indices12) = this%fqext (indices12) + this%nodes(i)%fqext
       
       call sparse_smatrix%assemblyIndexed(indices12 + qOffset, indices12 + qOffset, -0.5d0*this%nodes(i)%kextq)
       
    end do

    return
        
  end subroutine beamefk
  
end module class_beam
