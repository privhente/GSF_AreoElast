module class_shell

  use class_node_6
  use class_shell_element
  use class_sparse_matrix

  implicit none

  type :: shell
         
     integer :: nnodes
     integer :: nelements
     integer :: nncoordinates ! number of nodal coordinates = number of nodes times 6
     integer :: necoordinates ! number of enhanced coordinates = number of elements times 8
     integer :: tncoordinates ! total number of coordinates
     type(node_6), allocatable :: nodes(:)
     type(shell_element), allocatable :: elements(:)
     integer, allocatable :: indicesq6(:),indicesq8(:)
     real(kind = 8) :: ie, penergy, deltat
     integer :: simutype
     real(kind = 8), dimension(:), allocatable :: f, fqint, fqdyn, fv, fqext
     logical :: flag_kgeo_on, flag_kmat_on
     
   contains

     procedure :: initialization => shellinitialization
     procedure :: actualization => shellactualization
     procedure, private :: shellm
     procedure, private :: shellefk
     
  end type shell
  
  contains

subroutine  shellinitialization(this, sparse_mass, offset)

  implicit none
  
  class(shell), intent(inout) :: this
  class(sparse_matrix), intent(inout) :: sparse_mass
  integer, intent(in) :: offset

  integer :: i, i6, i8, j, j6

  this%nncoordinates = 6*this%nnodes

  this%necoordinates = 8*this%nelements

  this%tncoordinates = this%nncoordinates+this%necoordinates

  allocate(this%fqint (this%tncoordinates))
  allocate(this%fqdyn (this%tncoordinates))
  allocate(this%fv    (this%tncoordinates))
  allocate(this%fqext (this%tncoordinates))
  
  do i = 1, this%nnodes

     i6 = 6*(i-1)
     
     this%nodes(i)%coordinates(:) = [i6+1, i6+2, i6+3, i6+4, i6+5, i6+6]
     
  end do

  do i = 1, this%nelements
     
     i8 = 8*(i-1)+this%nncoordinates
     
     this%elements(i)%coordinates(:) = [i8+1, i8+2, i8+3, i8+4, i8+5, i8+6, i8+7, i8+8]
     
     do j = 1, 4 
        
        j6 = 6*(j-1)
        
        this%elements(i)%q_0(j6+1:j6+6) = this%nodes(this%elements(i)%connectivity(j))%q_0
        
     end do
     
     call this%elements(i)%precalculation()
     
  end do
  
  call this%shellm(sparse_mass, offset)

  return
  
end subroutine shellinitialization

subroutine shellactualization(this, q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2, v_1, v_2)

  implicit none
    
  class(shell), intent(inout) :: this
  real(kind = 8), intent(in) :: q_1(:)
  real(kind = 8), intent(in) :: ml_a(:)
  real(kind = 8), intent(in) :: sl_a(:)
  
  type(sparse_matrix),  intent(inout)  :: sparse_smatrix
  integer, intent(in) :: qOffset, vOffset
  
  real(kind = 8), optional, intent(in) :: q_2(:)
  real(kind = 8), optional, intent(in) :: v_1(:)
  real(kind = 8), optional, intent(in) :: v_2(:)
  
  if (present(q_2)) then

     if (present(v_1) .and. present(v_2)) then

        call this%shellefk(q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2, v_1, v_2)
        
     else
     
        call this%shellefk(q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2)

     end if
     
  else
     
     call this%shellefk(q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset)
     
  end if
  
  return
  
end subroutine shellactualization

!> @brief calculate mass matrix
!
!> generate mass matrix from elements in the shell
!> distribute terms according to element connectivity

subroutine shellm(this, sparse_mass, offset)
  implicit none
  
  class(shell), intent(inout) :: this 
  class(sparse_matrix), intent(inout) :: sparse_mass
  integer, intent(in) :: offset   
  
  integer :: i, indices24(24)
  
  do i = 1, this%nelements
     
     indices24 = [this%nodes(this%elements(i)%connectivity(1))%coordinates(:), &
                  this%nodes(this%elements(i)%connectivity(2))%coordinates(:), &
                  this%nodes(this%elements(i)%connectivity(3))%coordinates(:), &
                  this%nodes(this%elements(i)%connectivity(4))%coordinates(:)]

     call sparse_mass%assemblyIndexed(indices24 + offset, indices24 + offset, this%elements(i)%massmatrix())
     
  end do
    
  return
    
end subroutine shellm

!> @brief calculate nodal forces and stiffness matrix for shell object
!
!> updates forces and stiffness for current configuration and loading
!> writes resulting stiffness matrix directly to the iteration matrix
!
!> @param[in] q_1 Geometry configuration at beginning of time step (position and one director)
!> @param[in] ml_a material load vector at temporal midpoint (components 1:3 are forces, 4:6 are moments)
!> @param[in] sl_a spatial load vector at temporal midpoint (components 1:3 are forces, 4:6 are moments)
!> @param[in] smatrix iteration matrix
!> @param[in] qOffset offset to the position part of the iteration matrix
!> @param[in] vOffset offset to the velocity part of the iteration matrix
!> @param[in] q_2 Geometry configuration at current increment (position and director)
!> @param[in] v_1 Velocity configuration at beginning of time step (position and one director)
!> @param[in] v_2 Velocity configuration at current increment (position and director)

subroutine shellefk(this, q_1, ml_a, sl_a, sparse_smatrix, qOffset, vOffset, q_2, v_1, v_2)

    implicit none

    class(shell), intent(inout) :: this
    real(kind = 8), intent(in) :: q_1(:)
    real(kind = 8), intent(in) :: ml_a(:)
    real(kind = 8), intent(in) :: sl_a(:)
    
    type(sparse_matrix),  intent(inout)  :: sparse_smatrix
    integer, intent(in) :: qOffset, vOffset
    
    real(kind = 8), optional, intent(in) :: q_2(:)
    real(kind = 8), optional, intent(in) :: v_1(:)
    real(kind = 8), optional, intent(in) :: v_2(:)
    
    integer :: i, indices6(6), indices8(8), indices24(24), indices32(32)
        
    this%penergy   = 0.0d0
    this%fqint(:)  = 0.0d0
    this%fqdyn(:)  = 0.0d0
    this%fv(:)     = 0.0d0    
    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(indices24, indices8, indices32)
    
    do i = 1, this%nelements
       indices24 = [this%nodes(this%elements(i)%connectivity(1))%coordinates, &
                    this%nodes(this%elements(i)%connectivity(2))%coordinates, &
                    this%nodes(this%elements(i)%connectivity(3))%coordinates, &
                    this%nodes(this%elements(i)%connectivity(4))%coordinates]
       indices8 = [this%elements(i)%coordinates]
       indices32 = [indices24, indices8]

       this%elements(i)%deltat = this%deltat
       
       if (present(q_2)) then
          if (present(v_1) .and. present(v_2)) then
             call this%elements(i)%internalterms(q_1(indices32), q_2(indices32), v_1(indices32), v_2(indices32))
          else
             call this%elements(i)%internalterms(q_1(indices32), q_2(indices32))
          end if
       else
          call this%elements(i)%internalterms(q_1(indices32))
       end if
    end do
!$OMP  END PARALLEL DO
! this is formulated in terms of addition, so there will be data corruption if run in parallel
    do i = 1, this%nelements
       
       indices24 = [this%nodes(this%elements(i)%connectivity(1))%coordinates, &
                    this%nodes(this%elements(i)%connectivity(2))%coordinates, &
                    this%nodes(this%elements(i)%connectivity(3))%coordinates, &
                    this%nodes(this%elements(i)%connectivity(4))%coordinates]

       indices8 = [this%elements(i)%coordinates]
       
       this%penergy                   = this%penergy                  +this%elements(i)%penergy
       
       this%fqint (indices24)            = this%fqint (indices24)           +this%elements(i)%fqint ( 1:24)
       this%fqint (indices8 )            = this%fqint (indices8 )           +this%elements(i)%fqint (25:32)
       
       this%fqdyn (indices24)            = this%fqdyn (indices24)           +this%elements(i)%fqdyn ( 1:24)
       this%fqdyn (indices8 )            = this%fqdyn (indices8 )           +this%elements(i)%fqdyn(25:32)

       this%fv (indices24)            = this%fv (indices24)           +this%elements(i)%fv ( 1:24)
       this%fv (indices8 )            = this%fv (indices8 )           +this%elements(i)%fv (25:32)

       call sparse_smatrix%assemblyIndexed(indices24 + qOffset, indices24 + qOffset, this%elements(i)%kqq( 1:24,  1:24))
       call sparse_smatrix%assemblyIndexed(indices24 + qOffset, indices8  + qOffset, this%elements(i)%kqq( 1:24, 25:32))
       call sparse_smatrix%assemblyIndexed(indices8  + qOffset, indices24 + qOffset, this%elements(i)%kqq(25:32,  1:24))
       call sparse_smatrix%assemblyIndexed(indices8  + qOffset, indices8  + qOffset, this%elements(i)%kqq(25:32, 25:32))

       call sparse_smatrix%assemblyIndexed(indices24 + vOffset, indices24 + vOffset, this%elements(i)%kvv( 1:24,  1:24))
       call sparse_smatrix%assemblyIndexed(indices24 + vOffset, indices8  + vOffset, this%elements(i)%kvv( 1:24, 25:32))
       call sparse_smatrix%assemblyIndexed(indices8  + vOffset, indices24 + vOffset, this%elements(i)%kvv(25:32,  1:24))
       call sparse_smatrix%assemblyIndexed(indices8  + vOffset, indices8  + vOffset, this%elements(i)%kvv(25:32, 25:32))
       
       call sparse_smatrix%assemblyIndexed(indices24 + qOffset, indices24 + vOffset, this%elements(i)%kqv( 1:24,  1:24))
       call sparse_smatrix%assemblyIndexed(indices24 + qOffset, indices8  + vOffset, this%elements(i)%kqv( 1:24, 25:32))
       call sparse_smatrix%assemblyIndexed(indices8  + qOffset, indices24 + vOffset, this%elements(i)%kqv(25:32,  1:24))
       call sparse_smatrix%assemblyIndexed(indices8  + qOffset, indices8  + vOffset, this%elements(i)%kqv(25:32, 25:32))
       
       call sparse_smatrix%assemblyIndexed(indices24 + vOffset, indices24 + qOffset, this%elements(i)%kvq( 1:24,  1:24))
       call sparse_smatrix%assemblyIndexed(indices24 + vOffset, indices8  + qOffset, this%elements(i)%kvq( 1:24, 25:32))
       call sparse_smatrix%assemblyIndexed(indices8  + vOffset, indices24 + qOffset, this%elements(i)%kvq(25:32,  1:24))
       call sparse_smatrix%assemblyIndexed(indices8  + vOffset, indices8  + qOffset, this%elements(i)%kvq(25:32, 25:32))
              
    end do
    

        
! evaluation of external terms at the nodes, material and spatial loads.

    this%fqext (:)    = 0.0d0
    
    do i = 1, this%nnodes
       
       indices6 = this%nodes(i)%coordinates

       if (present(q_2)) then
          
          call this%nodes(i)%externalterms(0.5d0*(q_1(indices6)+q_2(indices6)), ml_a(indices6), sl_a(indices6))
          
       else

          call this%nodes(i)%externalterms(       q_1(indices6)               , ml_a(indices6), sl_a(indices6))
          
       end if
       
       this%fqext (indices6) = this%fqext (indices6) + this%nodes(i)%fqext
       
       call sparse_smatrix%assemblyIndexed(indices6 + qOffset, indices6 + qOffset, -0.5d0*this%nodes(i)%kextq)
              
    end do
    
    return
        
  end subroutine shellefk
  
end module class_shell
