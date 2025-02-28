module class_model_structure
  use class_node_6
  use class_node_12
  use class_rigid_body
  use class_beam
  use class_shell
  use class_constraint_6
  use class_constraint_12
  use class_constraint_6_to_12
  use class_boundary
  use class_load_6
  use class_load_12
  use class_point_mass
  use class_sparse_matrix
  use class_amplitude
  use class_spring12
  use class_damping12
  use class_nullspace
  use my_math_structure, only: eye, vec21mat36symm, vec6mat9symm
  use class_matrix12

  implicit none

  type :: step
    character(len = :), allocatable :: simutype
    real(kind = 8)                  :: totalt, deltat, tolerance
    real(kind = 8)                  :: emin, emax
    integer                         :: iterationlimit
    integer                         :: gravityflag
    integer                         :: arc_len_method_flag, number_desired_iteration
    integer                         :: nEF
    integer                         :: step
    integer                         :: i_simutype
    integer                         :: sparse_flag, output_flag, optional_adaptive_deltat
  end type step

  type :: post_proc
    integer, allocatable :: id(:)
    character(:), allocatable :: strtype
    logical :: boolean_e = .TRUE., boolean_q = .TRUE., boolean_v = .TRUE., boolean_m = .TRUE., boolean_lambda = .TRUE.
  end type post_proc

  type :: model_structure
     ! global node lists
     type(node_6), allocatable :: nodes6(:) !< nodes with 6 coordinates for shells
     type(node_12), allocatable :: nodes12(:) !< nodes with 12 coordinates for rigid bodies and beams

     ! structures
     type(rigid_body), allocatable :: bodies(:) !< rigid bodies
     type(beam), allocatable :: beams(:) !< beams
     type(shell), allocatable :: shells(:) !< shells

     ! constraints
     type(constraint_6_c), allocatable :: constraints6(:) !< constraints involving nodes with 6 coordinates
     type(constraint_12_c), allocatable :: constraints12(:) !< constraints involving nodes with 12 coordinates
     type(constraint_6to12_c), allocatable :: constraints6to12(:) !< constraints involving nodes6 and node12!!! with 12 coordinates
     
     ! boundaries
     type(amplitude), allocatable :: boundaryamplitudes6(:) !< amplitudes involving constraints with 6 coordinates
     type(amplitude), allocatable :: boundaryamplitudes12(:) !< amplitudes involving constraints with 12 coordinates
     type(boundary_12), allocatable :: boundaries12(:) !< boundary transformation for constraints with 12 coordinates
     type(boundary_6),  allocatable ::  boundaries6(:) !< boundary transformation for constraints with  6 coordinates

     ! condensed interactions and springs
     type(spring12_c), allocatable :: spring12(:)   !< springs for nodes with 12 coordinates
     type(damping12_c), allocatable :: damping12(:)  !< damping element for nodes with 12 coordinates
     
     ! coupling element
     integer :: nnodes_matrix12 !< total number of matrices for nodes with 12 coordinates
     integer :: number_matrix12 !< total number of different matrices 
     type(matrix_12) :: matrix12 !< coupling element for nodes with 12 coordinates
     
     ! loads
     type(amplitude), allocatable :: loadamplitudes6(:) !< load amplitudes involving nodes with 6 coordinates
     type(amplitude), allocatable :: loadamplitudes12(:) !< load amplitudes involving nodes with 12 coordinates
     type(load_6), allocatable :: loads6(:) !< loads involving nodes with 6 coordinates
     type(load_12), allocatable :: loads12(:) !< loads involving nodes with 6 coordinates

     type(step), allocatable :: settings(:)
     type(post_proc), allocatable :: postproc(:)

     ! point mass properties
     type(point_mass), allocatable :: pmass6(:) !< point mass for nodes with 6 coordinates.
     type(point_mass), allocatable :: pmass12(:) !< point mass for nodes with 12 coordinates.
     integer :: nnodes6 !< number of nodes with 6 coordinates
     integer :: nnodes12 !< number of nodes with 12 coordinates
     integer :: nbodies !< number of rigid bodies
     integer :: nbeams !< number of beams
     integer :: nshells !< number of shells
     integer :: nconstraints6 !< number of constraints for nodes with 6 coordinates.
     integer :: nconstraints12 !< number of constraints for nodes with 12 coordinates.
     integer :: nconstraints6to12 !< number of constraints for nodes with 12 coordinates.
     integer :: nboundaries12 !< number of constraints for nodes with 6 coordinates.
     integer :: nboundaries6  !< number of boundary conditions with 6 coordinates.
     integer :: nloads6 !< number of loads for nodes with 6 coordinates.
     integer :: nloads12 !< number of loads for nodes with 12 coordinates.
     integer :: npmass6 !< number of point masses for nodes with 6 coordinates.
     integer :: npmass12 !< number of point masses for nodes with 12 coordinates.
     integer :: nspring12 !< number of springs for nodes with 12 coordinates.
     integer :: nspring12properties !< number of springs properties for nodes with 12 coordinates.
     integer :: ndamping12 !< number of springs for nodes with 12 coordinates.
     integer :: ndamping12properties !< number of springs properties for nodes with 12 coordinates.
     integer :: trconstraints6 !< total rank of the constraints due to node6.
     integer :: trconstraints12 !< total rank of the constraints due to node12.
     integer :: trconstraints6to12 !< total rank of the constraints due to node12.
     integer :: trconstraints !< total rank of the constraints, including node6 and node12.
     integer :: tncoordinates6 !< total number of coordinates due to node6
     integer :: tncoordinates8 !< total number of coordinates due to internal degrees of freedom for shells.
     integer :: tncoordinates12 !< total number of coordinates due to node12.
     integer :: tncoordinates !< total number of coordinates, including node6 and node12 and strain dofs.
     integer :: tnstress    !< total number of elements in the system, exluded rigid bodies
     integer :: tnnullspace !< total number of null space if the constraint matrix
     integer :: tnnodes !< total number of nodes, including node6 and node12.
     type(sparse_matrix) :: sparse_mass !< the (constant) sparse mass matrix.
     type(nullspace) :: null_space !< null space of constraint matrix

     real(kind = 8), allocatable :: g(:)          !< constraints.
     real(kind = 8), allocatable :: q_0(:)        !< initial coordinates
     real(kind = 8), allocatable :: q_t(:)        !< coordinates at the end of analysis
     real(kind = 8), allocatable :: qdot_0(:)     !< initial time derivative of the coordinates
     real(kind = 8), allocatable :: qdot_t(:)     !< velocity
     real(kind = 8), allocatable :: lambda_t(:)   !< Lagrange multiplier.
     real(kind = 8), allocatable :: lambda_0(:)  !< Lagrange multiplier.
     
     real(kind = 8) :: lmomentum_t(3), amomentum_t(3)
     real(kind = 8) :: kenergy_t, strainenergy_t, gravenergy_t, penergy_t, tenergy_t
     real(kind = 8) :: kenergy, strainenergy, gravenergy, penergy, tenergy
     real(kind = 8) :: lmomentum(3), amomentum(3)

     real(kind = 8) :: time_t
     real(kind = 8) :: deltat

     real(kind = 8), allocatable :: fqint(:) !< internal forces associated with q.
     real(kind = 8), allocatable :: fqdyn(:) !< dynamic forces associated with velocity.
     real(kind = 8), allocatable :: fqext(:) !< external forces.
     real(kind = 8), allocatable :: fqext0(:) !< prescribed external forces vector, i.e. for aerodynamical loads
     real(kind = 8), allocatable :: fv(:) !< forces associated with velocities.
     real(kind = 8), allocatable :: fqgrav(:) !< forces associated with gravity.
     real(kind = 8), allocatable :: fqL(:) !< forces associated with constraint forces to coordinates.

     real(kind = 8), allocatable :: fqint_t(:) !< internal forces associated with q.
     real(kind = 8), allocatable :: fqdyn_t(:) !< dynamic forces associated with velocity.
     real(kind = 8), allocatable :: fqext_t(:) !< external forces.
     real(kind = 8), allocatable :: fqext0_t(:) !< external forces.
     real(kind = 8), allocatable :: fv_t(:) !< forces associated with velocities.
     real(kind = 8), allocatable :: fqL_t(:) !< forces associated with constraint forces.
     real(kind = 8) :: gravity(3) !< acceleration vector due to gravity
     real(kind = 8), allocatable :: stress(:),stress_t(:)
     logical :: boolean_abort = .FALSE.
     logical :: boolean_model_read = .FALSE.
     logical :: boolean_oldsim = .FALSE.
     logical :: boolean_write_init_state = .TRUE.
     logical :: boolean_outputfiles_open = .FALSE.
     character(:), allocatable :: output_filename, oldsimfilename

   contains
     procedure :: localSimuSetting
     procedure :: constructor
     procedure :: modelfk
     procedure :: modelgdg
     procedure :: modelndgq
     procedure :: modelloads
     procedure :: modelboundaries
     procedure :: modelfgrav
     procedure :: modelfL
     procedure :: readinginputs
     procedure :: read_oldsim
     procedure :: modelinvariants
     procedure :: output_open
     procedure :: output_write
     procedure :: output_close
     procedure :: output_write_model
     procedure, private :: modelm
     procedure, private :: constraintsrankcounter
     procedure, private :: constraintlist
  end type model_structure

contains

  subroutine constructor(this)
    implicit none

    class(model_structure), intent(inout) :: this

    integer :: i, i6, i12, j, j6, j8, j12, k, incoordinates, incoordinates6, incoordinates8, incoordinates12, instress, inodes_gl
    integer :: ilast, inull, aa, rank, rank_i
    integer :: node_i, n_constraints, n_nodes12, n_nodes6
    integer :: r1
    logical :: boolean_node
    integer :: c12, c6
    integer :: npoints
    real(kind = 8) :: slope
    logical :: boolean_nullspace
    
 

! allocating and initializing
    print*, 'Allocating and initializing elements and bodies...'
    incoordinates6  = 0
    incoordinates8  = 0
    incoordinates12 = 0
    instress        = 0
    
! rigid bodies
    if (allocated(this%bodies)) then
       do i = 1, this%nbodies
          incoordinates12 = incoordinates12+12
       end do
    end if

! beams
    if (allocated(this%beams)) then
       do i = 1, this%nbeams
          incoordinates12 = incoordinates12+12*this%beams(i)%nnodes
          instress        = instress + 6*this%beams(i)%nelements
       end do
    end if

! shells
    if (allocated(this%shells)) then
       do i = 1, this%nshells
          incoordinates6 = incoordinates6+6*this%shells(i)%nnodes
          incoordinates8 = incoordinates8+8*this%shells(i)%nelements
          !instress       = instress + 6*this%shells(i)%nelements
       end do
    end if

! coordinates
    this%tncoordinates6  = incoordinates6
    this%tncoordinates8  = incoordinates8
    this%tncoordinates12 = incoordinates12
    this%tncoordinates   = incoordinates6+incoordinates8+incoordinates12
    this%tnstress        = instress
    
    if (this%tncoordinates .eq. 0) then
      this%boolean_abort = .True.
      return
    end if
    
! nodes
    print*, 'Allocating nodes...'
    this%nnodes12 = 0
    this%nnodes6 = 0

    if (allocated(this%bodies)) then
       do i = 1, this%nbodies
          this%nnodes12 = this%nnodes12+1
       end do
    end if

    if (allocated(this%beams)) then
       do i = 1, this%nbeams
          this%nnodes12 = this%nnodes12+this%beams(i)%nnodes
       end do
    end if

    if (allocated(this%shells)) then
       do i = 1, this%nshells
          this%nnodes6 = this%nnodes6+this%shells(i)%nnodes
       end do
    end if
    this%tnnodes = this%nnodes6+this%nnodes12

! constraints to get rank of constraints
    call this%constraintsrankcounter()

! allocating more things
    if (this%nnodes6 > 0) then
       allocate(this%nodes6 (this%nnodes6))
    end if

    if (this%nnodes12 > 0) then
       allocate(this%nodes12(this%nnodes12))
    end if

    allocate(this%q_0    (this%tncoordinates))
    allocate(this%q_t    (this%tncoordinates))
    allocate(this%qdot_0 (this%tncoordinates))
    allocate(this%qdot_t (this%tncoordinates))
    allocate(this%lambda_0 (this%trconstraints))
    allocate(this%lambda_t (this%trconstraints))

    this%sparse_mass = sparse_matrix(this%tncoordinates, this%tncoordinates)
    allocate(this%g(this%trconstraints))
    allocate(this%fqint  (this%tncoordinates))
    allocate(this%fqdyn  (this%tncoordinates))
    allocate(this%fqext  (this%tncoordinates))
    allocate(this%fqext0 (this%tncoordinates))
    allocate(this%fqL    (this%tncoordinates))
    allocate(this%stress (this%tnstress))

    allocate(this%fv     (this%tncoordinates))
    allocate(this%fqgrav (this%tncoordinates))

    allocate(this%fqint_t(this%tncoordinates))
    allocate(this%fqdyn_t(this%tncoordinates))
    allocate(this%fqext_t(this%tncoordinates))
    allocate(this%fqext0_t(this%tncoordinates))
    allocate(this%fqL_t  (this%tncoordinates))
    allocate(this%fv_t   (this%tncoordinates))
    allocate(this%stress_t(this%tnstress))

! assigning initial conditions, in... = initial 
    incoordinates = 0
    instress      = 0
    inodes_gl     = 0
    
! rigid bodies
    if (allocated(this%bodies)) then
       do i = 1, this%nbodies
          this%bodies(i)%indicesq12(1:12)               = (/ (r1, r1=12*(i-1)+1, 12*(i-1)+12)  /)
          this%bodies(i)%node%node_structure%globalID   = i
          this%q_0   (incoordinates+1:incoordinates+12) = this%bodies(i)%node%q_0
          this%q_t   (incoordinates+1:incoordinates+12) = this%bodies(i)%node%q_0
          this%qdot_0(incoordinates+1:incoordinates+12) = 0.0d0
          this%qdot_t(incoordinates+1:incoordinates+12) = 0.0d0
          incoordinates = incoordinates+12
          inodes_gl     = i
       end do
    end if

! beams
    if (allocated(this%beams)) then
       do i = 1, this%nbeams
          allocate(this%beams(i)%indicesq12(12*this%beams(i)%nnodes))
          do j = 1, this%beams(i)%nnodes
             j12 = 12*(j-1)
	           this%beams(i)%indicesq12(j12+1:j12+12) = (/ (r1, r1=incoordinates+j12+1,incoordinates+j12+12)  /)
             this%q_0   (incoordinates+j12+1:incoordinates+j12+12) = this%beams(i)%nodes(j)%q_0
             this%q_t   (incoordinates+j12+1:incoordinates+j12+12) = this%beams(i)%nodes(j)%q_0
             this%qdot_0(incoordinates+j12+1:incoordinates+j12+12) = 0.0d0
             this%qdot_t(incoordinates+j12+1:incoordinates+j12+12) = 0.0d0
             
             this%beams(i)%nodes(j)%node_structure%globalID = inodes_gl + j
          end do
          
          !< global indices for beam element stress vector
          allocate(this%beams(i)%indicesst(6*this%beams(i)%nelements))
          this%beams(i)%indicesst(1:6*this%beams(i)%nelements) = (/ (r1, r1=instress+1,instress+6*this%beams(i)%nelements)  /)

          incoordinates = incoordinates + 12*this%beams(i)%nnodes
          inodes_gl     = inodes_gl + this%beams(i)%nnodes
          instress      = instress + 6*this%beams(i)%nelements
       end do
    end if

! shells
    if (allocated(this%shells)) then
       do i = 1, this%nshells
          allocate(this%shells(i)%indicesq6(6*this%shells(i)%nnodes))
          do j = 1, this%shells(i)%nnodes ! nodal coordinates.
             j6 = 6*(j-1)
             this%shells(i)%indicesq6(j6+1:j6+6) = (/ (r1, r1=incoordinates+j6+1,incoordinates+j6+6)  /)
             this%q_0   (incoordinates+j6+1:incoordinates+j6+6) = this%shells(i)%nodes(j)%q_0 ! nodal
             this%q_t   (incoordinates+j6+1:incoordinates+j6+6) = this%shells(i)%nodes(j)%q_0 ! nodal
             this%qdot_0(incoordinates+j6+1:incoordinates+j6+6) = 0.0d0 ! nodal
             this%qdot_t(incoordinates+j6+1:incoordinates+j6+6) = 0.0d0 ! nodal
             
            this%shells(i)%nodes(j)%node_structure%globalID = inodes_gl + j             
          end do
          incoordinates = incoordinates+6*this%shells(i)%nnodes
          inodes_gl     = inodes_gl + this%shells(i)%nnodes
          
          allocate(this%shells(i)%indicesq8(this%shells(i)%nelements*8))
          do j = 1, this%shells(i)%nelements ! elemental coordinates.
             j8 = 8*(j-1)
             r1 = 0
             this%shells(i)%indicesq8(j8+1:j8+8) = (/ (r1, r1=incoordinates+j8+1,incoordinates+j8+8)  /)
             this%q_0   (incoordinates+j8+1:incoordinates+j8+8) = 0.0d0 ! elemental
             this%q_t   (incoordinates+j8+1:incoordinates+j8+8) = 0.0d0 ! elemental
             this%qdot_0(incoordinates+j8+1:incoordinates+j8+8) = 0.0d0 ! elemental
             this%qdot_t(incoordinates+j8+1:incoordinates+j8+8) = 0.0d0 ! elemental
          end do
          incoordinates = incoordinates+8*this%shells(i)%nelements
       end do
    end if

 ! Allocating constant stiffness and load ranges for spring12
    if (allocated(this%spring12)) then
      do i = 1,this%nspring12
        npoints = size(this%spring12(i)%c%data_points,1)
        allocate(this%spring12(i)%c%ranges(npoints+1,2))
        allocate(this%spring12(i)%c%const_stiff(npoints))

        this%spring12(i)%c%ranges      = 0.0d0
        this%spring12(i)%c%const_stiff = 0.0d0

        do j = 1,npoints-1
          this%spring12(i)%c%ranges(j+1,1)    = (this%spring12(i)%c%data_points(j,1)+this%spring12(i)%c%data_points(j+1,1))/2
          slope = (this%spring12(i)%c%data_points(j+1,2)-this%spring12(i)%c%data_points(j,2))/(this%spring12(i)%c%data_points(j+1,1)-this%spring12(i)%c%data_points(j,1))
          this%spring12(i)%c%ranges(j+1,2)    = slope*(this%spring12(i)%c%ranges(j+1,1)-this%spring12(i)%c%data_points(j,1))+this%spring12(i)%c%data_points(j,2)
          this%spring12(i)%c%const_stiff(j)   = slope
        end do
        this%spring12(i)%c%ranges(npoints+1,1)  = this%spring12(i)%c%data_points(npoints,1)
        this%spring12(i)%c%ranges(npoints+1,2)  = this%spring12(i)%c%data_points(npoints,2)
        this%spring12(i)%c%const_stiff(npoints) = slope
      end do
  end if
 !
  ! Allocating constant stiffness and load ranges for spring12
    if (allocated(this%damping12)) then
      do i = 1,this%ndamping12
        npoints = size(this%damping12(i)%c%data_points,1)
        allocate(this%damping12(i)%c%ranges(npoints+1,2))
        allocate(this%damping12(i)%c%const_stiff(npoints))

        this%damping12(i)%c%ranges      = 0.0d0
        this%damping12(i)%c%const_stiff = 0.0d0

        do j = 1,npoints-1
          this%damping12(i)%c%ranges(j+1,1)    = (this%damping12(i)%c%data_points(j,1)+this%damping12(i)%c%data_points(j+1,1))/2
          slope = (this%damping12(i)%c%data_points(j+1,2)-this%damping12(i)%c%data_points(j,2))/(this%damping12(i)%c%data_points(j+1,1)-this%damping12(i)%c%data_points(j,1))
          this%damping12(i)%c%ranges(j+1,2)    = slope*(this%damping12(i)%c%ranges(j+1,1)-this%damping12(i)%c%data_points(j,1))+this%damping12(i)%c%data_points(j,2)
          this%damping12(i)%c%const_stiff(j)   = slope
        end do
        this%damping12(i)%c%ranges(npoints+1,1)  = this%damping12(i)%c%data_points(npoints,1)
        this%damping12(i)%c%ranges(npoints+1,2)  = this%damping12(i)%c%data_points(npoints,2)
        this%damping12(i)%c%const_stiff(npoints) = slope
      end do
    end if

! Assigning global indizes for nodes 12
  ilast = 0
  if (allocated(this%nodes12)) then
       do i = 1, this%nnodes12
          i12 = 12*(i-1)
          this%nodes12(i)%coordinates = [i12+ 1, i12+ 2, i12+ 3, i12+ 4, i12+ 5, i12+ 6, i12+ 7, i12+ 8, i12+ 9, i12+10, i12+11, i12+12]      ! indices for coordinates
          this%nodes12(i)%velocity    = [i12+ 1, i12+ 2, i12+ 3, i12+ 4, i12+ 5, i12+ 6, i12+ 7, i12+ 8, i12+ 9, i12+10, i12+11, i12+12]      ! indices for velocity, HERE: same node dofs for the velocity!
       end do
       ilast = 12*this%nnodes12
    end if

! Assigning global indizes for nodes 6
    if (allocated(this%nodes6)) then
       k  = 0
       i6 = ilast
       do i = 1, this%nshells
          do j = 1, this%shells(i)%nnodes
             j6 = 6*(j-1)+i6
             this%nodes6(k+j)%coordinates   = [j6+1, j6+2, j6+3, j6+4, j6+5, j6+6]
          end do
          k = k+this%shells(i)%nnodes
          ! This line is for over-jumping the strain degrees of freedom
          i6 = i6+6*this%shells(i)%nnodes+8*this%shells(i)%nelements
       end do
       ilast = j6+6

       ! assigning global indizes for node 6 but without strain dof
       k = 0
       do i = 1, this%nshells
          do j = 1, this%shells(i)%nnodes
            k = k + 1
	          this%nodes6(k)%coordinates_q(1:6) = 12*this%nnodes12 + (/ (r1, r1 = 6*(k-1)+1,6*(k-1)+6)  /)
          end do
       end do

    end if

! Check, if nullspace is required for simulation
    boolean_nullspace = .FALSE.
    do i = 1, size(this%settings)
        if (this%settings(i)%simuType .eq. 'modal' .or. this%settings(i)%simuType .eq. 'buckling') then
          boolean_nullspace = .TRUE.
        end if
    end do

! Creating a constraint list for null-space calculation: This is used to perform a modal analysis or linear buckling analysis
   this%tnnullspace = 0
   if (boolean_nullspace .eqv. .TRUE.) then
      call this%constraintlist()      ! Rearanging the constraint lists, with merging common constraints
      inull  = 0
      rank   = 0
      if (allocated(this%null_space%list)) then
        do i = 1,size(this%null_space%list,1)

            ! Determine number of rank and indices for merged constraints
            rank_i = 0
            n_constraints = 0
            if (allocated(this%null_space%list(i)%constraints)) then
              n_constraints = size(this%null_space%list(i)%constraints)
            end if
            
            do j = 1,n_constraints

              c12 = count(this%null_space%constraint_list(this%null_space%list(i)%constraints(j),1:2) .GT. 0)
              c6  = count(this%null_space%constraint_list(this%null_space%list(i)%constraints(j),3:4) .GT. 0)

              ! checking, which kind of constraint
              if (c12 > 0 .AND. c6 > 0)         then
                rank_i = rank_i + this%constraints6to12(this%null_space%constraint_list(this%null_space%list(i)%constraints(j),5))%c%rankcount
              else if (c12 > 0 .AND. c6 == 0) then
                rank_i = rank_i + this%constraints12(this%null_space%constraint_list(this%null_space%list(i)%constraints(j),5))%c%rankcount
              elseif (c12 == 0 .AND. c6 > 0) then
               rank_i = rank_i + this%constraints6(this%null_space%constraint_list(this%null_space%list(i)%constraints(j),5))%c%rankcount
              end if

            end do

            this%null_space%list(i)%rank = rank_i
            
            this%null_space%list(i)%n_nodes12 = 0
            if (allocated(this%null_space%list(i)%nodes12)) then
              this%null_space%list(i)%n_nodes12 = size(this%null_space%list(i)%nodes12)
            end if
            n_nodes12 = this%null_space%list(i)%n_nodes12
            
            this%null_space%list(i)%n_nodes6 = 0
            if (allocated(this%null_space%list(i)%nodes6)) then
              this%null_space%list(i)%n_nodes6 = size(this%null_space%list(i)%nodes6)
            end if
            n_nodes6 = this%null_space%list(i)%n_nodes6
            
            ! Loop to get local indices for the nodes
            allocate(this%null_space%list(i)%node(n_nodes12 + n_nodes6))
            aa = 0
            do j = 1,n_nodes12
              aa = aa + 1
              allocate(this%null_space%list(i)%node(aa)%indices_temp_q(12))
              this%null_space%list(i)%node(aa)%indices_temp_q(1:12) = (/ (r1, r1 = 12*(j-1)+1,12*(j-1)+12)  /)
            end do

            do j = 1,n_nodes6
              aa = aa + 1
              allocate(this%null_space%list(i)%node(aa)%indices_temp_q(6))
              this%null_space%list(i)%node(aa)%indices_temp_q(1:6) = 12*n_nodes12 + (/ (r1, r1 = 6*(j-1)+1,6*(j-1)+6) /)
            end do

            ! Loop to get null-space indices of node 12 and node 6
            if (rank_i > 12*n_nodes12 + 6*n_nodes6) then
              print *,'Error: Check constraint set',this%null_space%list(i)%constraints(:)
              !stop
            elseif (rank_i < 12*n_nodes12 + 6*n_nodes6)  then
              allocate(this%null_space%list(i)%indicesn(12*n_nodes12 + 6*n_nodes6 - rank_i))
              this%null_space%list(i)%rank = rank_i
              do j = 1, (12*n_nodes12 + 6*n_nodes6 - rank_i)
                inull = inull + 1
                this%null_space%list(i)%indicesn(j) = inull
              end do
            end if
            rank = rank + rank_i
        end do
      end if
      
      ! Extra loop to consider the rest nodes6 in null space matrix
      aa = 0
      if (allocated(this%nodes6)) then
        do i = 1,size(this%nodes6)
          node_i = i
          boolean_node = .FALSE. !.TRUE.
          do j = 1,size(this%null_space%list,1)
              if (ANY(this%null_space%list(j)%nodes6 .eq. node_i)) then
                boolean_node = .TRUE.
                exit
              else
                boolean_node = .FALSE.
              end if
          end do
          if (boolean_node .eqv. .FALSE.) then
            aa = aa + 1
          end if
        end do

        if (aa .ne. 0) then
          allocate(this%null_space%addlist6(aa))    ! addlist6 = [node_i, indices_q(nod_i), indicesn]
          allocate(this%null_space%ns_add6(aa))
          aa = 0
          do i = 1,size(this%nodes6)
            node_i = i
            boolean_node = .FALSE.
            do j = 1,size(this%null_space%list,1)
              if (ANY(this%null_space%list(j)%nodes6 .eq. node_i)) then
                boolean_node = .TRUE.
                exit
              else
                boolean_node = .FALSE.
              end if
            end do
            if (boolean_node .eqv. .FALSE.) then
                aa = aa + 1
                this%null_space%addlist6(aa) = node_i
                allocate(this%null_space%ns_add6(aa)%indicesn(6))
	              this%null_space%ns_add6(aa)%indicesn(1:6) = (/ (r1, r1 = inull+1,inull+6) /)
                inull = inull + 6
            end if
          end do
        end if
      end if
      ! Total number of null spaces
      this%tnnullspace = this%tncoordinates6 + this%tncoordinates12 - rank
    end if
    

    ! Defining global Mass matrix
    call this%modelm()

    ! Initialize global vectors and initial setting
    call this%localSimuSetting(1.0d0, 0,.TRUE.,.TRUE.)
    this%time_t      = 0.0d0
    this%lambda_t    = 0.0d0
    this%q_t         = this%q_0
    this%qdot_t      = this%qdot_0
 
    this%fqint_t  =  0.0d0
    this%fqdyn_t  =  0.0d0
    this%fqext_t  =  0.0d0
    this%fqext0_t =  0.0d0
    this%fqL_t    =  0.0d0
    this%fv_t     =  0.0d0
    this%fqext0   =  0.0d0
    this%stress   =  0.0d0
    
    this%lmomentum_t    = 0.0d0
    this%amomentum_t    = 0.0d0
    this%kenergy_t      = 0.0d0
    this%strainenergy_t = 0.0d0
    this%gravenergy_t   = 0.0d0
    this%penergy_t      = 0.0d0
    this%tenergy_t      = 0.0d0
    this%stress_t       = 0.0d0
    
    this%lmomentum    = 0.0d0
    this%amomentum    = 0.0d0
    this%kenergy      = 0.0d0
    this%strainenergy = 0.0d0
    this%gravenergy   = 0.0d0
    this%penergy      = 0.0d0
    this%tenergy      = 0.0d0
    
    return
  end subroutine constructor

!< subroutine for calculating energies and linear/angular momentum
  subroutine modelinvariants(this, q1, v1_opt)
    implicit none
    class(model_structure), intent(inout) :: this
    real(kind = 8), intent(in) :: q1(:)
    real(kind = 8), intent(in), optional :: v1_opt(:)
    real(kind = 8), allocatable :: p1(:), v1(:)
    real(kind = 8) :: amomentum_m(3), d_amomentum(3), rotation(3,3)
    integer :: i, j, k, i12, i6, j6
    
    allocate(p1(this%tncoordinates))
    allocate(v1(this%tncoordinates))
    
    v1 = 0.0d0
    if (present(v1_opt)) v1 = v1_opt
    
    !< global momentum
    call mkl_dcsrmv('N', this%tncoordinates, this%tncoordinates, 1.0d0, 'G  F  ', &
         this%sparse_mass%csr_values, this%sparse_mass%csr_columns, &
         this%sparse_mass%csr_rowIndices(1:size(this%sparse_mass%csr_rowIndices) - 1), &
         this%sparse_mass%csr_rowIndices(2:size(this%sparse_mass%csr_rowIndices)), &
         v1, 0.0d0, p1)

    !< energy balance
    this%kenergy      = 0.5d0*dot_product(v1, p1)
    this%strainenergy = this%penergy
    this%gravenergy   = dot_product(this%fqgrav,q1-this%q_0)
    this%penergy      = -(this%gravenergy - this%strainenergy) !< Potential energy = strain energy + energy for gravitational forces
    this%tenergy      = this%kenergy + this%penergy
       
    !< calcuation of momentum     
    this%lmomentum(:) = 0.0d0
    this%amomentum(:) = 0.0d0
    amomentum_m(:)    = 0.0d0
    d_amomentum(:)    = 0.0d0
    rotation(:,:)     = 0.0d0
    
    do i = 1, this%nbodies
        i12 = 12*(i-1)
        this%lmomentum  = this%lmomentum + p1(i12+1:i12+3)
        d_amomentum     = cross(q1(i12+1:i12+3), p1(i12+1:i12+3)) + cross(q1(i12+4:i12+6), p1(i12+4:i12+6)) + cross(q1(i12+7:i12+9), p1(i12+7:i12+9)) + cross(q1(i12+10:i12+12), p1(i12+10:i12+12))            
        this%amomentum  = this%amomentum + d_amomentum
        rotation(:,:)   = outer(q1(i12+4:i12+6),this%q_0(i12+4:i12+6)) + outer(q1(i12+7:i12+9),this%q_0(i12+7:i12+9)) + outer(q1(i12+10:i12+12),this%q_0(i12+10:i12+12))
        amomentum_m     = amomentum_m + matmul(transpose(rotation),d_amomentum)
    end do
    
    do i = this%nbodies+1, this%nnodes12
        i12 = 12*(i-1)
        this%lmomentum = this%lmomentum + p1(i12+ 1:i12+3)
        d_amomentum    = cross(q1(i12+1:i12+3), p1(i12+1:i12+3)) + &
                          cross(q1(i12+4:i12+6), p1(i12+4:i12+6)) + &
                          cross(q1(i12+7:i12+9), p1(i12+7:i12+9))
        this%amomentum = this%amomentum + d_amomentum 
        rotation(:,:)  = outer(q1(i12+4:i12+6),this%q_0(i12+4:i12+6)) + outer(q1(i12+7:i12+9),this%q_0(i12+7:i12+9)) + outer(q1(i12+10:i12+12),this%q_0(i12+10:i12+12))
        amomentum_m    = amomentum_m + matmul(transpose(rotation),d_amomentum)
    end do
    
    k = 0
    i6 = 12*this%nnodes12
    do i = 1, this%nshells
        do j = 1, this%shells(i)%nnodes
          j6 = 6*(j-1)+i6         
          this%lmomentum = this%lmomentum + p1(j6+ 1:j6+ 3)
          d_amomentum    = cross(q1(j6+ 1:j6+ 3), p1(j6+ 1:j6+ 3)) + cross(q1(j6+ 4:j6+ 6), p1(j6+ 4:j6+ 6))
          this%amomentum = this%amomentum + d_amomentum
        end do
        k = k+this%shells(i)%nnodes
        i6 = i6+6*this%shells(i)%nnodes+8*this%shells(i)%nelements
    end do

    return
  end subroutine modelinvariants
  
!< Subroutine for opening result files
  subroutine output_open(this, flag_opt)
    implicit none

    class(model_structure), intent(inout) :: this
    integer(8), parameter :: block_size = 2 ** 200
    integer :: i, j, io_error
    character(:), allocatable :: output_filename
    character(len=7) :: str_status = 'replace'
    character(len=6) :: str_position = 'rewind'
    integer, optional, intent(in) :: flag_opt
    integer :: flag
    
3232 format(1000000E30.15e3)
3233 format(1000000i10)

    flag = 0
    if (present(flag_opt)) flag = flag_opt

    output_filename = trim(adjustl(this%output_filename))

    if (this%boolean_oldsim) then
      str_status   = 'unknown'
      str_position = 'append'
    end if
    
    if (allocated(this%oldsimfilename)) then
      if (this%output_filename .eq. this%oldsimfilename) this%boolean_write_init_state = .FALSE.
    end if

    !< open output files
    if (flag .eq. 1) then 
      open(unit = 2000, file = output_filename // '_steps.dres', action = 'readwrite', status = str_status, POSITION = str_position)
      if (this%boolean_write_init_state) write(2000,3233) this%settings(1)%step, this%settings(1)%i_simutype, 1, 0
    end if
    if (flag .eq. 2) then
      open(unit = 2100, file = output_filename // '_t.dres', action = 'readwrite', status = str_status, POSITION = str_position)
      if (this%boolean_write_init_state) write(2100,3232) this%time_t, 0.0d0
    end if
    if (flag .eq. 3) then
      open(unit = 3000, file = output_filename // '_q.dres', action = 'readwrite', status = str_status, POSITION = str_position) !, BUFFERED='YES', BLOCKSIZE=block_size)
      if (this%boolean_write_init_state) write(3000,3232) this%q_t(:)
    end if
    if (flag .eq. 4) then
      open(unit = 3100, file = output_filename // '_v.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      if (this%boolean_write_init_state) write(3100,3232) this%qdot_t(:)
    end if
    if (flag .eq. 5) then 
      open(unit = 3200, file = output_filename // '_lambda.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      if (this%boolean_write_init_state) write(3200,3232) this%lambda_t
    end if
    if (flag .eq. 6) then
      open(unit = 3300, file = output_filename // '_e.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      if (this%boolean_write_init_state) write(3300,3232) this%kenergy_t, this%penergy_t, this%tenergy_t
    end if
    if (flag .eq. 7) then
      open(unit = 3400, file = output_filename // '_m.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      if (this%boolean_write_init_state) write(3400,3232) this%lmomentum_t(:), this%amomentum_t(:)
    end if
    
    if (flag .eq. 8) then
      open(unit = 3500, file = output_filename // '_stress.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      if (this%boolean_write_init_state) write(3500,3232) this%stress_t(:)
    end if
    
    !if (flag .eq. 9) then
    !  open(unit = 3500, file = output_filename // '_ce.txt', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    !  if (this%boolean_write_init_state) write(3600,3232) this%stress_t(:)
    !end if
        
    if (flag .eq. 0) then
      open(unit = 2000, file = output_filename // '_steps.dres', action = 'readwrite', status = str_status, POSITION = str_position)
      open(unit = 2100, file = output_filename // '_t.dres', action = 'readwrite', status = str_status, POSITION = str_position)
      open(unit = 3000, file = output_filename // '_q.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      open(unit = 3100, file = output_filename // '_v.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      open(unit = 3200, file = output_filename // '_lambda.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      open(unit = 3300, file = output_filename // '_e.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      open(unit = 3400, file = output_filename // '_m.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      open(unit = 3500, file = output_filename // '_stress.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      !open(unit = 3600, file = output_filename // '_ce.dres', action = 'readwrite', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
      
      if (this%boolean_write_init_state) then
        write(2000,3233) this%settings(1)%step, this%settings(1)%i_simutype, 1, 0
        write(2100,3232) this%time_t, 0.0d0
        write(3000,3232) this%q_t(:)
        write(3100,3232) this%qdot_t(:)
        write(3200,3232) this%lambda_t
        write(3300,3232) this%kenergy_t, this%penergy_t, this%tenergy_t
        write(3400,3232) this%lmomentum_t(:), this%amomentum_t(:)
        write(3500,3232) this%stress_t(:)
        !write(3600,3232) this%time_t, this%matrix12%matrix12_object(1)%lie_u(:), this%matrix12%matrix12_object(1)%num(:), this%matrix12%matrix12_object(1)%ua_dot(:)
      end if
      this%boolean_outputfiles_open = .TRUE.
    end if
    
  end subroutine output_open

!< subroutine for writing model to file
  subroutine output_write_model(this)
    implicit none
    class(model_structure), intent(in) :: this
    integer :: i, j, io_error
    character(:), allocatable :: output_filename
    
    output_filename = trim(adjustl(this%output_filename))
    
    !< write model_structure structure
    open(unit = 999, file = output_filename // '_model.m', action = 'readwrite', status = 'replace')
      if (allocated(this%bodies)) then
        write(999, *) 'nbodies = ', size(this%bodies), ';'
      end if

      if (allocated(this%beams)) then
        write(999, *) 'nbeams = ', size(this%beams), ';'
        do i = 1, size(this%beams)
          write(999, *) 'model.beams(', i, ').nnodes = '   , this%beams(i)%nnodes   , ';'
          write(999, *) 'model.beams(', i, ').nelements = ', this%beams(i)%nelements, ';'
          do j = 1, this%beams(i)%nelements
            write(999, *) 'model.beams(', i, ').connectivity(', j, ', :) = [', this%beams(i)%elements(j)%connectivity(:), '];'
          end do
        end do
      end if

      if (allocated(this%shells)) then
        write(999, *) 'nshells = ', size(this%shells), ';'
        do i = 1, size(this%shells)
          write(999, *) 'model.shells(', i, ').nnodes = ', this%shells(i)%nnodes, ';'
          write(999, *) 'model.shells(', i, ').nelements = ', this%shells(i)%nelements, ';'
          do j = 1, this%shells(i)%nelements
            write(999, *) 'model.shells(', i, ').connectivity(', j, ', :) = ['
            write(999, *) this%shells(i)%elements(j)%connectivity(:), '];'
          end do
        end do
      end if
    close(999)
    
    return
  end subroutine output_write_model
  
!< Subroutine for writing result files
  subroutine output_write(this, step, t, subtime, istep, i_simuType, iteration, q, v, lambda, stress, flag_opt)
    implicit none
    class(model_structure), intent(in) :: this
    real(kind = 8), intent(in) :: t, subtime
    real(kind = 8), dimension(:), intent(in) :: q, v, lambda, stress
    integer, intent(in) :: step, i_simuType, istep, iteration
    integer, optional, intent(in) :: flag_opt
    integer :: flag
    
    flag = 0
    if (present(flag_opt)) flag = flag_opt
  
3232 format(1000000E30.15e3)
3233 format(1000000i10)
    
    if (flag .eq. 1) write(2000,3233) step, i_simuType, istep, iteration
    if (flag .eq. 2) write(2100,3232) t, subtime
    if (flag .eq. 3) write(3000,3232) q(:)
    if (flag .eq. 4) write(3100,3232) v(:)
    if (flag .eq. 5) write(3200,3232) lambda(:)
    if (flag .eq. 6) write(3300,3232) this%kenergy, this%penergy, this%tenergy
    if (flag .eq. 7) write(3400,3232) this%lmomentum(:), this%amomentum(:)
    if (flag .eq. 8) write(3500,3232) this%stress(:)
    !if (flag .eq. 9) write(3600,3232) t, this%matrix12%matrix12_object(1)%lie_u(:), this%matrix12%matrix12_object(1)%num(:), this%matrix12%matrix12_object(1)%ua_dot(:)

    if (flag .eq. 0) then
      write(2000,3233) step, i_simuType, istep, iteration
      write(2100,3232) t, subtime
      write(3000,3232) q(:)
      write(3100,3232) v(:)
      write(3200,3232) lambda(:)
      write(3300,3232) this%kenergy, this%penergy, this%tenergy
      write(3400,3232) this%lmomentum(:), this%amomentum(:)
      write(3500,3232) this%stress(:)
      !write(3600,3232)  t, this%matrix12%matrix12_object(1)%lie_u(:), this%matrix12%matrix12_object(1)%num(:), this%matrix12%matrix12_object(1)%ua_dot(:)
    end if
    
    return
  end subroutine output_write

  ! Subroutine for writing result files
  subroutine output_close(this, flag_opt)
    implicit none
    class(model_structure), intent(in) :: this
    integer, optional, intent(in) :: flag_opt
    integer :: flag
    
    flag = 0
    if (present(flag_opt)) flag = flag_opt
    if (flag .eq. 1) close(2000)
    if (flag .eq. 2) close(2100)
    if (flag .eq. 3) close(3000)
    if (flag .eq. 4) close(3100)
    if (flag .eq. 5) close(3200)
    if (flag .eq. 6) close(3300)
    if (flag .eq. 7) close(3400)
    if (flag .eq. 8) close(3500)
    !if (flag .eq. 9) close(3600)
    
    if (flag .eq. 0) then
      close(2000)
      close(2100)
      close(3000)
      close(3100)
      close(3200)
      close(3300)
      close(3400)
      close(3500)
      !close(3600)
    end if
    
    return
  end subroutine output_close
  
!< subroutine for local simulation settings for shell and beam elements  
  subroutine localSimuSetting(this, deltat, simutype,flag_kgeo_on,flag_kmat_on)
    implicit none
    class(model_structure), intent(inout) :: this
    real(kind = 8), intent(in)  :: deltat
    integer, intent(in) :: simutype
    logical, intent(in) :: flag_kgeo_on, flag_kmat_on

    integer :: i, j
    
    if (allocated(this%bodies)) then
      do i = 1, this%nbodies
        this%bodies(i)%deltat = deltat
        this%bodies(i)%simutype = simutype
      end do
    end if
    
    if (allocated(this%beams)) then
      do i = 1, this%nbeams
        this%beams(i)%deltat  = deltat
        this%beams(i)%simutype = simutype
        this%beams(i)%flag_kgeo_on = flag_kgeo_on
        this%beams(i)%flag_kmat_on = flag_kmat_on
        do j = 1, size(this%beams(i)%elements)
          this%beams(i)%elements(j)%deltat = deltat
          this%beams(i)%elements(j)%simutype = simutype
          this%beams(i)%elements(j)%flag_kgeo_on = flag_kgeo_on
          this%beams(i)%elements(j)%flag_kmat_on = flag_kmat_on
        end do
      end do
    end if
    
    if (allocated(this%shells)) then
      do i = 1, this%nshells
        this%shells(i)%deltat = deltat
        this%shells(i)%simutype = simutype
        this%shells(i)%flag_kgeo_on = flag_kgeo_on
        this%shells(i)%flag_kmat_on = flag_kmat_on
        do j = 1, size(this%shells(i)%elements)
          this%shells(i)%elements(j)%deltat = deltat
          this%shells(i)%elements(j)%simutype = simutype
          this%shells(i)%elements(j)%flag_kgeo_on = flag_kgeo_on
          this%shells(i)%elements(j)%flag_kmat_on = flag_kmat_on
        end do
      end do
    end if
    
    return
  end subroutine
  
subroutine constraintlist(this)
    implicit none
    class(model_structure), intent(inout) :: this

    integer, allocatable :: list(:,:), list_nodes12(:,:), list_nodes6(:,:)
    integer, allocatable :: flagList(:)
    integer, allocatable :: constraint_list(:,:)
    integer, allocatable :: vecNodes12(:), vecNodes6(:)
    integer :: flag12, flag6
    integer :: a, a12, a6, j, i, b, c, k, bmax

  ! initializing: Merging all constraint into a single list, sorted by node 12 and node 6 constraints
    allocate(constraint_list(this%nconstraints6to12 + this%nconstraints12 + this%nconstraints6,5))
    allocate(flagList(this%nconstraints6to12 + this%nconstraints12 + this%nconstraints6))

    flagList = 0

    i = 0
    constraint_list (:,:) = 0
    if (allocated(this%constraints6to12)) then
      do while (i /= this%nconstraints6to12)
        i = i + 1
        constraint_list(i,1) = this%constraints6to12(i)%c%constraint%nodes(1) ! node12
        constraint_list(i,2) = 0
        constraint_list(i,3) = this%constraints6to12(i)%c%constraint%nodes(2) ! node6
        constraint_list(i,4) = 0
        constraint_list(i,5) = i
      end do
    end if

    j = 0
    if (allocated(this%constraints12)) then
      do while (j /= this%nconstraints12)
        i = i + 1
        j = j + 1
        constraint_list(i,1) = this%constraints12(j)%c%constraint%nodes(1) ! node12 - node 1
        constraint_list(i,2) = this%constraints12(j)%c%constraint%nodes(2) ! node12 - node 2
        constraint_list(i,3) = 0
        constraint_list(i,4) = 0
        constraint_list(i,5) = j
      end do
    end if
    
    j = 0
    if (allocated(this%constraints6)) then
      do while (j /= this%nconstraints6)
        i = i + 1
        j = j + 1
        constraint_list(i,1) = 0
        constraint_list(i,2) = 0
        constraint_list(i,3) = this%constraints6(j)%c%constraint%nodes(1) ! node6 - node 1
        constraint_list(i,4) = this%constraints6(j)%c%constraint%nodes(2) ! node6 - node 2
        constraint_list(i,5) = j
      end do
    end if
    
  ! Loop over all constraints to find dependent constraints for null-space calculation
    j    = 0
    a    = 0
    bmax = 0
    do while (j /= size(constraint_list,1))
      j = j + 1
      i = 0
      b = 0
      if (flagList(j) == 0) then
        a = a + 1
        ! find entries unequal zero
        c = count(constraint_list(j,1:2) .GT. 0)
        allocate(vecnodes12(c))
        if (c /= 0) then
          vecnodes12 = constraint_list(j,1:c)
        end if

        c = count(constraint_list(j,3:4) .GT. 0)
        allocate(vecnodes6(c))
        if (c /= 0) then
          vecnodes6 = constraint_list(j,3:2+c)
        end if

        ! Loop to merge constraints according to equal nodes, result (=list) is an index-list for constraint_list
        do while (i /= size(constraint_list,1))
          i = i + 1
          if (flagList(i) == 0) then
            flag12 = 0
            if (size(vecNodes12) .ne. 0) then
              call search_add_List(constraint_list(i,1:2), vecNodes12, flag12)  ! search in constraint_list for all nodes in vecNodes12
              if (flag12 == 1) then
                if (constraint_list(i,3) /= 0) then
                  call append_int_vec(vecNodes6, constraint_list(i,3)) ! append node nr to vecNodes6, in case it isn't inside
                end if
              end if
            end if
            flag6 = 0
            if (size(vecNodes6) .ne. 0) then
              call search_add_List(constraint_list(i,3:4), vecNodes6,  flag6)
              if (flag6 == 1) then
                if (constraint_list(i,1) /= 0) then
                  call append_int_vec(vecNodes12, constraint_list(i,1))  ! append node nr to vecNodes12, in case it isn't inside
                end if
              end if
            end if
            if (flag12 == 1 .or. flag6 == 1) then
              b = b + 1
              if (b > bmax) then
                bmax = b
              end if
              call reallocate_IArray(list,a,b,i)
              flagList(i) = 1
              i = 0
            end if
          end if
        end do

        do k = 1,size(vecnodes12)
            if (k==1) then
              call reallocate_IArray(list_nodes12,a,1,a)
            end if
            call reallocate_IArray(list_nodes12,a,k+1,vecnodes12(k))
        end do
        deallocate(vecnodes12)

        do k = 1,size(vecnodes6)
            if (k==1) then
              call reallocate_IArray(list_nodes6,a,1,a)
            end if
            call reallocate_IArray(list_nodes6,a,k+1,vecnodes6(k))
        end do
        deallocate(vecnodes6)

      end if

    end do

    ! Loop to fill the global constraint list, which is needed for the null_space calculation
    allocate(this%null_space%list(a))
    do i = 1, a
      allocate(this%null_space%list(i)%constraints(count(list(i,:) .GT. 0)))
      this%null_space%list(i)%constraints = list(i,1:count(list(i,:) .GT. 0))

      if (allocated(list_nodes12)) then
        if (i .le. size(list_nodes12,1)) then
          if (list_nodes12(i,1) .eq. i) then
              allocate(this%null_space%list(i)%nodes12(count(list_nodes12(i,:) .GT. 0)-1))
              this%null_space%list(i)%nodes12 = list_nodes12(i,2:count(list_nodes12(i,:) .GT. 0))
          end if
        end if
      else
        allocate(this%null_space%list(i)%nodes12(0))
      end if

      if (allocated(list_nodes6)) then
        if (list_nodes6(i,1) .eq. i) then
          allocate(this%null_space%list(i)%nodes6(count(list_nodes6(i,:) .GT. 0)-1))
          this%null_space%list(i)%nodes6 = list_nodes6(i,2:count(list_nodes6(i,:) .GT. 0))
        end if
      else
        allocate(this%null_space%list(i)%nodes6(0))
      end if
    end do

    allocate(this%null_space%constraint_list (size(constraint_list,1),size(constraint_list,2)))
    this%null_space%constraint_list = constraint_list
    
    return
end subroutine constraintlist

! subroutine to search same nodes and add nodes to the list
subroutine search_add_List(vec1, vecNodes, flag)
    implicit none

    integer :: flag, i, j, val1, val2
    integer, allocatable, intent(inout) :: vecNodes(:)
    integer, intent(in) :: vec1(:)
    integer, allocatable :: vec2(:), temp_nodes(:)

    flag = 0
    allocate(vec2(size(vecNodes)))
    vec2 = vecNodes
    ! First search, if search node is in current (constraint) list
    do i = 1, size(vec1)
      val1 = vec1(i)
      do j = 1, size(vec2)
        val2 = vec2(j)
        if (val1 == val2) then
          flag = 1
          exit
        end if
      end do
      if (flag == 1) then
        exit
      end if
    end do

    if (flag == 1) then
      allocate(temp_nodes(size(vec1) + size(vec2)))
      temp_nodes(1:size(vec1))                         = vec1
      temp_nodes(size(vec1)+1:size(vec1) + size(vec2)) = vecNodes
      do i = 1,size(temp_nodes)
        if (temp_nodes(i) /= 0) then
            if (.NOT. ANY(vecNodes .eq. temp_nodes(i))) then
              call append_int_vec(vecNodes,temp_nodes(i))
            end if
        end if
      end do
    end if

end subroutine search_add_List

!< Subroutine to append data to array, but only, if the value isn't already inside
subroutine append_int_vec(vec, val)
    implicit none
    integer, allocatable, intent(inout) :: vec(:)
    integer, intent(in) :: val
    integer, allocatable :: temp_vec(:)
    integer :: num_add

    if (.not. allocated(vec)) then
      allocate(vec(1))
      vec(1) = -1
      num_add = 0
    else
      num_add = 1
    end if

    if (ALL(vec .ne. val)) then
      allocate(temp_vec(size(vec) + num_add))
      temp_vec(1:size(vec)) = vec
      temp_vec(size(vec) + num_add) = val
      deallocate(vec)
      allocate(vec(size(temp_vec)))
      vec = temp_vec
    end if

  end subroutine

! Subroutine to reallocate a vector with type integer
subroutine reallocate_ivec(vec_a,vec_b)
    implicit none
    integer, allocatable, intent(inout) :: vec_a(:)
    integer, allocatable :: temp_vec_a(:)
    integer :: vec_b(:)
    integer :: n_a, n_b

    n_a      = size(vec_a)
    if (.not. allocated(vec_a)) then
      n_a = 0
    end if

    n_b = size(vec_b)
    if (n_a /= 0) then
      allocate(temp_vec_a(n_a))
      temp_vec_a = vec_a
      deallocate(vec_a)
      allocate(vec_a(n_a+n_b))
      vec_a(1:n_a) = temp_vec_a(:)
    else
      allocate(vec_a(n_a+n_b))
    end if
    vec_a(n_a+1:n_a+n_b) = vec_b(:)

  end subroutine reallocate_ivec

subroutine modelm(this)
    implicit none

    class(model_structure), intent(inout) :: this
    integer :: incoordinates, liq, i, i6, i12
    integer :: indices12(12), indices3(3)
    integer :: indices_q_total(24)  ! Indices of the two nodes coupled with the coupling element 

    indices12 = (/1,2,3,4,5,6,7,8,9,10,11,12 /)
    indices3 = (/1,2,3/)
    incoordinates = 0

! rigid bodies
    if (allocated(this%bodies)) then
       print*, 'Initializing rigid bodies...'
       do i = 1, this%nbodies
          call this%bodies(i)%initialization() ! Was passiert hier? Wo ist die subroutine initialization? 
          liq = incoordinates
          call this%sparse_mass%assemblyIndexed(indices12 + liq, indices12 + liq, this%bodies(i)%m)
          incoordinates = incoordinates+12
       end do
    end if

! coupling element
    if (allocated(this%matrix12%matrix12_object)) then
        print*, 'Initializing coupling elements'
        do i = 1, this%nnodes_matrix12
            indices_q_total(13:24) = this%nodes12(this%matrix12%matrix12_object(i)%node2)%coordinates
           
            !< checking if mass matrix is located between two nodes or connected to the environment
            if (this%matrix12%matrix12_object(i)%node1 .eq. 0) then 
                indices_q_total(1:12) = indices_q_total(13:24)
            else 
                indices_q_total(1:12) = this%nodes12(this%matrix12%matrix12_object(i)%node1)%coordinates
            end if 
            
            !< building the mass matrix
            call  this%matrix12%matrix12_built_mass(i)
             
            !< assamble the added mass matrix into the total sparse matrix 
             if (this%matrix12%matrix12_object(i)%node1 .eq. 0) then 
                 call this%sparse_mass%assemblyIndexed(indices_q_total(1:12), indices_q_total(1:12), this%matrix12%matrix12_object(i)%M_tilde(1:12,1:12))
             else 
                 call  this%sparse_mass%assemblyIndexed(indices_q_total, indices_q_total, this%matrix12%matrix12_object(i)%M_tilde)
             end if 
             
        end do 
    end if
        
! beams
    if (allocated(this%beams)) then
       print*, 'Initializing beams...'
       do i = 1, this%nbeams
          liq = incoordinates
          call this%beams(i)%initialization(this%sparse_mass, liq)
          incoordinates = incoordinates+this%beams(i)%ncoordinates
       end do
    end if

! shells
    if (allocated(this%shells)) then
       print*, 'Initializing shells...'
       do i = 1, this%nshells
          liq = incoordinates
          call this%shells(i)%initialization(this%sparse_mass, liq)
          incoordinates = incoordinates+this%shells(i)%nncoordinates+this%shells(i)%necoordinates
       end do
    end if

!! pointmass12
    if (allocated(this%pmass12)) then
       print*, 'Initializing point masses 12...'
       do i = 1, this%npmass12
          i12 = 12*(this%pmass12(i)%node-1)
          call this%sparse_mass%assemblyIndexed(i12 + indices3, i12 + indices3, this%pmass12(i)%mass * eye(3))
       end do
    end if

!! pointmass6
    if (allocated(this%pmass6)) then
       print*, 'Initializing point masses 6...'
       do i = 1, this%npmass6
          i6 = 6*(this%pmass6(i)%node-1)+12*this%nnodes12
          call this%sparse_mass%assemblyIndexed(i6 + indices3, i6 + indices3, this%pmass6(i)%mass * eye(3))
       end do
    end if
    
    ! data is now in sparse matrix, initialize
    if (this%sparse_mass%nrows .ne. 0) then
      print*, 'Calculating sparse mass matrix...'
      call this%sparse_mass%initSparseStorage()
    end if 

  return
end subroutine modelm

subroutine modelfL(this,lambda,sparse_smatrix)
    implicit none

    class(model_structure), intent(inout) :: this
    type(sparse_matrix),  intent(inout)   :: sparse_smatrix

    real(kind = 8) :: lambda(:)

    integer :: row, index, col
    integer :: tncoordinates

    tncoordinates = this%tncoordinates

  ! Calculate constraint forces. Perform calculation dm0gm0_t * lambda. This is index-wise the matrix Kqg * lambda
    this%fqL(:) = 0.0d0
    do row = 1, tncoordinates
      do index = sparse_smatrix%csr_rowIndices(row), sparse_smatrix%csr_rowIndices(row + 1) - 1
        col = sparse_smatrix%csr_columns(index)
        if (col > 2*tncoordinates) then
          this%fqL(row) = this%fqL(row) + lambda(col - 2*tncoordinates) * sparse_smatrix%csr_values(index)
        endif
      end do
    end do

  end subroutine modelfL

subroutine modelfgrav(this)

    implicit none

    class(model_structure), intent(inout) :: this

    integer :: j12, i, j, j6
    integer :: incoordinates, tncoordinates

    real(kind = 8), dimension(:), allocatable :: gravity_vector

    tncoordinates = this%tncoordinates
    
    if (allocated(gravity_vector) .eqv. .FALSE.) allocate(gravity_vector(tncoordinates))

! computing gravity forces if they are present.
    print *, "Initializing gravity..."

    incoordinates = 0

    gravity_vector(:) = 0.0d0

! rigid bodies

    if (allocated(this%bodies)) then

      do i = 1, this%nbodies

          gravity_vector(incoordinates+1:incoordinates+3) = this%gravity

          incoordinates = incoordinates+12

      end do

    end if

! beams

    if (allocated(this%beams)) then

      do i = 1, this%nbeams

          do j = 1, this%beams(i)%nnodes

            j12 = 12*(j-1)

            gravity_vector(incoordinates+j12+1:incoordinates+j12+3) = this%gravity

          end do

          incoordinates = incoordinates+12*this%beams(i)%nnodes

      end do

    end if

! shells

    if (allocated(this%shells)) then

      do i = 1, this%nshells

          do j = 1, this%shells(i)%nnodes

            j6 = 6*(j-1)

            gravity_vector(incoordinates+j6+1:incoordinates+j6+3) = this%gravity

          end do

          incoordinates = incoordinates+6*this%shells(i)%nnodes+8*this%shells(i)%nelements

      end do

    end if

    call mkl_dcsrmv('N', tncoordinates, tncoordinates, 1.0d0, 'G  F  ', &
      this%sparse_mass%csr_values, this%sparse_mass%csr_columns, &
      this%sparse_mass%csr_rowIndices(1:size(this%sparse_mass%csr_rowIndices) - 1), &
      this%sparse_mass%csr_rowIndices(2:size(this%sparse_mass%csr_rowIndices)), &
      gravity_vector, 0.0d0, this%fqgrav)

  end subroutine modelfgrav

!> @brief update model_structure forces and stiffness
!
!> Assembles stiffness block matrices (qq, qv, vq, vv) directly into the given iteration matrix
!
!> @param[in] q_1 Geometry configuration at beginning of timestep
!> @param[in] q_2 Geometry configuration at current iteration
!> @param[in] v_1 Velocities at beginning of timestep
!> @param[in] v_2 Velocities at current iteration
!> @param[inout] smatrix The iteration matrix
!> @param[in] indicesv Velocity indices for assembly
!> @param[in] indicesq Positional indices for assembly
!> @param[in] ml_t Material loading at current time
!> @param[in] sl_t Spatial force loading at current time
subroutine modelfk(this, q_1, ml_t, sl_t, sparse_smatrix, indicesq, indicesv, q_2, v_1, v_2)  
    implicit none

    class(model_structure), intent(inout) :: this
    real(kind = 8), intent(in) :: q_1(:)
    real(kind = 8), intent(in) :: ml_t(:)
    real(kind = 8), intent(in) :: sl_t(:)

    type(sparse_matrix),  intent(inout)  :: sparse_smatrix
    integer, allocatable, dimension(:),  intent(in) :: indicesq, indicesv

    real(kind = 8), optional, intent(in) :: q_2(:)
    real(kind = 8), optional, intent(in) :: v_1(:)
    real(kind = 8), optional, intent(in) :: v_2(:)

    integer :: innodes6, innodes12, incoordinates, liq, uiq, i, lil, uil
    integer :: indices12a_q(12), indices12b_q(12), indices12a_v(12), indices12b_v(12), indices12a_q_b(12),  indices12a_v_a(12), indices12a_v_b(12)
    integer :: indices_q_total(24), indices12v_q_total(24) !Indices of the two nodes coupled by the coupling element    
    real(kind = 8), dimension(21) :: va0
    
    ! HACK: dummy parameter to make PARDISO work
    real(kind = 8), allocatable :: dummy_mat(:,:)

    allocate(dummy_mat(1,1))
    dummy_mat = tiny(0.0d0)

    this%penergy = 0.0d0
    this%fqint   = 0.0d0
    this%fqdyn   = 0.0d0
    this%fqL     = 0.0d0
    this%fv      = 0.0d0
    this%fqext   = 0.0d0
    this%stress  = 0.0d0
    
    innodes6 = 0
    innodes12 = 0
    incoordinates = 0
    
    ! HACK: zeroes to the main diagonal to make PARDISO work
    do i = 1, size(indicesq)
      call sparse_smatrix%assemblyIndexed(indicesq(i:i), indicesq(i:i), dummy_mat)
      call sparse_smatrix%assemblyIndexed(indicesv(i:i), indicesv(i:i), dummy_mat)
    end do

! rigid bodies
    if (allocated(this%bodies)) then
       do i = 1, size(this%bodies)
          liq = incoordinates+1
          uiq = incoordinates+12
          lil = 6*innodes12+1
          uil = 6*(innodes12+1)
          
          this%bodies(i)%deltat = this%deltat
          
          if (present(q_2)) then
             if (present(v_1) .and. present(v_2)) then
                call this%bodies(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), q_2(liq:uiq), v_1(liq:uiq), v_2(liq:uiq))
             else
                call this%bodies(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), q_2(liq:uiq))
             end if
          else
             call this%bodies(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil))
          end if

          this%penergy         = this%penergy+this%bodies(i)%penergy
          this%fqint(liq:uiq)  = this%bodies(i)%fqint
          this%fqdyn(liq:uiq)  = this%bodies(i)%fqdyn
          this%fqext(liq:uiq)  = this%bodies(i)%fqext
          this%fv(liq:uiq)     = this%bodies(i)%fv

!          ! write stiffness block matrices to the appropriate indices in the iteration matrix
          call sparse_smatrix%assemblyIndexed(indicesq(liq:uiq), indicesq(liq:uiq), this%bodies(i)%kqq)
          call sparse_smatrix%assemblyIndexed(indicesv(liq:uiq), indicesv(liq:uiq), this%bodies(i)%kvv)
          call sparse_smatrix%assemblyIndexed(indicesq(liq:uiq), indicesv(liq:uiq), this%bodies(i)%kqv)
          call sparse_smatrix%assemblyIndexed(indicesv(liq:uiq), indicesq(liq:uiq), this%bodies(i)%kvq)
          innodes12     = innodes12    +1
          incoordinates = uiq
       end do
    end if

! beams
    if (allocated(this%beams)) then
       do i = 1, size(this%beams)
          liq = incoordinates+1
          uiq = incoordinates+this%beams(i)%ncoordinates
          
          lil = 6*innodes12+1
          uil = 6*(innodes12 +this%beams(i)%nnodes)
          
          this%beams(i)%deltat = this%deltat
          
          if (present(q_2)) then
             if (present(v_1) .and. present(v_2)) then
                call this%beams(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), sparse_smatrix, indicesq(liq)-1, indicesv(liq)-1, q_2(liq:uiq), v_1(liq:uiq), v_2(liq:uiq))
             else
                call this%beams(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), sparse_smatrix, indicesq(liq)-1, indicesv(liq)-1, q_2(liq:uiq))
             end if
          else
             call this%beams(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), sparse_smatrix, indicesq(liq)-1, indicesv(liq)-1)
          end if

          this%penergy               = this%penergy+this%beams(i)%penergy
          this%fqint(liq:uiq)        = this%beams(i)%fqint
          this%fqdyn(liq:uiq)        = this%beams(i)%fqdyn
          this%fqext(liq:uiq)        = this%beams(i)%fqext
          this%fv(liq:uiq)           = this%beams(i)%fv
       
          this%stress(this%beams(i)%indicesst) = this%beams(i)%stress
 
          innodes12     = innodes12 + this%beams(i)%nnodes
          incoordinates = uiq
          
       end do
    end if

! shells
    if (allocated(this%shells)) then
       do i = 1, size(this%shells)
          liq = incoordinates+1
          uiq = incoordinates+this%shells(i)%nncoordinates+this%shells(i)%necoordinates
          lil = 6*(innodes12+innodes6)+1
          uil = 6*(innodes12+innodes6+this%shells(i)%nnodes)
          
          this%shells(i)%deltat = this%deltat
          
          if (present(q_2)) then
             if (present(v_1) .and. present(v_2)) then
                call this%shells(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), sparse_smatrix, indicesq(liq)-1, indicesv(liq)-1, q_2(liq:uiq), v_1(liq:uiq), v_2(liq:uiq))
             else
                call this%shells(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), sparse_smatrix, indicesq(liq)-1, indicesv(liq)-1, q_2(liq:uiq))
             end if
          else
             call this%shells(i)%actualization(q_1(liq:uiq), ml_t(lil:uil), sl_t(lil:uil), sparse_smatrix, indicesq(liq)-1, indicesv(liq)-1)
          end if

          this%penergy        = this%penergy+this%shells(i)%penergy
          this%fqint(liq:uiq) = this%shells(i)%fqint
          this%fqdyn(liq:uiq) = this%shells(i)%fqdyn
          this%fqext(liq:uiq) = this%shells(i)%fqext
          this%fv(liq:uiq)    = this%shells(i)%fv

          innodes6      = innodes6 + this%shells(i)%nnodes
          incoordinates = uiq
       end do
    end if

! spring12
    if (allocated(this%spring12)) then
      do i = 1, this%nspring12
        indices12a_q = this%nodes12(this%spring12(i)%c%node)%coordinates
        
        this%spring12(i)%c%deltat = this%deltat
        
        call this%spring12(i)%c%spring12_f(q_2(indices12a_q), q_1(indices12a_q), this%q_0(indices12a_q))
        this%fqint(indicesq(indices12a_q)) = this%fqint(indices12a_q) + this%spring12(i)%c%fqint

        call this%spring12(i)%c%spring12_k(q_2(indices12a_q), q_1(indices12a_q), this%q_0(indices12a_q))
        call sparse_smatrix%assemblyIndexed(indicesq(indices12a_q), indicesq(indices12a_q), this%spring12(i)%c%kqq)

        this%penergy = this%penergy + this%spring12(i)%c%penergy
      end do
    end if

! coupling element  
    if (allocated(this%matrix12%matrix12_object)) then
        
        do i = 1, this%nnodes_matrix12
        this%matrix12%matrix12_object(i)%deltat = this%deltat
        
        indices12a_q_b = this%nodes12(this%matrix12%matrix12_object(i)%node2)%coordinates  
        indices12a_v_b = this%nodes12(this%matrix12%matrix12_object(i)%node2)%velocity 
               
        ! If node1 (A) = 0, there are no indices that can be assigned to node(0) and they are not needed  
        if (this%matrix12%matrix12_object(i)%node1 .ne. 0) then
            indices12a_q = this%nodes12(this%matrix12%matrix12_object(i)%node1)%coordinates !Indices of the node. Continuous for all nodes. 1st node: 1:12, 2nd: 13:24, ...
            indices12a_v_a = this%nodes12(this%matrix12%matrix12_object(i)%node1)%velocity 
        end if
        
        select case (this%matrix12%matrix12_object(i)%matrix12_strain_flag) 
            
        case('linear') 
                
                !actualization of the elastic forces 
                call this%matrix12%matrix12_lin_stiffness(q_2(indices12a_q_b), q_1(indices12a_q_b), this%q_0(indices12a_q_b), this%number_matrix12, i) 
                this%fqint(indicesq(indices12a_q_b)) = this%fqint(indicesq(indices12a_q_b)) + this%matrix12%matrix12_object(i)%fqint
                call sparse_smatrix%assemblyIndexed(indicesq(indices12a_q_b), indicesq(indices12a_q_b), this%matrix12%matrix12_object(i)%kqq)
        
                this%penergy = this%penergy + this%matrix12%matrix12_object(i)%penergy
        
        case('objective')  
               
                    !< First node = Second node at t = 0
                if (this%matrix12%matrix12_object(i)%node1 .eq. 0) then  
                    va0 = 0.0d0
            
                    !< Calculation of the elastic forces
                    call this%matrix12%matrix12_obj_stiffness(this%q_0(indices12a_q_b), this%q_0(indices12a_q_b), this%q_0(indices12a_q_b), q_2(indices12a_q_b), q_1(indices12a_q_b), this%q_0(indices12a_q_b),va0, va0, v_2(indices12a_v_b), v_1(indices12a_v_b), this%number_matrix12, i)
                    !< Calculation of the dynamic forces
                    call this%matrix12%matrix12_built_mass(i)
                    call this%matrix12%matrix12_mass(this%deltat, this%q_0(indices12a_q_b), this%q_0(indices12a_q_b), this%q_0(indices12a_q_b), q_2(indices12a_q_b), q_1(indices12a_q_b), this%q_0(indices12a_q_b), va0, va0, v_2(indices12a_v_b), v_1(indices12a_v_b), this%number_matrix12, i)
            
                    !< Actualization of the forces
                    !< Actualization of the elastic forces
                    this%fqint(indices12a_q_b) = this%fqint(indices12a_q_b) + this%matrix12%matrix12_object(i)%fqint_obj(13:24)
                    !< Actualization of dynamic forces
                    this%fqdyn(indices12a_q_b) = this%fqdyn(indices12a_q_b) + this%matrix12%matrix12_object(i)%fqdyn_AB(13:24) 
                    !< Actualization of the forces associated with velocities
                    this%fv(indices12a_v_b)    = this%fv(indices12a_v_b) + this%matrix12%matrix12_object(i)%fv_AB(13:24) 
          
                    !< Actualization of the iteration matrix
                    !< kqq (resulting from stiffnes terms) 
                    call sparse_smatrix%assemblyIndexed(indices12a_q_b, indices12a_q_b, this%matrix12%matrix12_object(i)%kqq_obj(13:24,13:24)) 
                    !< kvv (resulting from inertia terms) 
                    call sparse_smatrix%assemblyIndexed(indicesv(indices12a_v_b), indicesv(indices12a_v_b), this%matrix12%matrix12_object(i)%kvv(13:24,13:24))
                    !< kqv (resulting from inertia terms) 
                    call sparse_smatrix%assemblyIndexed(indices12a_q_b, indicesv(indices12a_v_b), this%matrix12%matrix12_object(i)%kqv(13:24,13:24))
                    !< kvq (resulting from inertia terms) 
                    call sparse_smatrix%assemblyIndexed(indicesv(indices12a_v_b), indices12a_q_b, this%matrix12%matrix12_object(i)%kvq(13:24,13:24))

                    !< Actualization of the total energy
                    this%penergy = this%penergy + this%matrix12%matrix12_object(i)%penergy   
                    
                else     
                    
                    indices_q_total(1:12)     = indices12a_q
                    indices_q_total(13:24)    = indices12a_q_b
                    indices12v_q_total(1:12)  = indices12a_v_a
                    indices12v_q_total(13:24) = indices12a_v_b
                    
                    !< Calling subroutine to calculate the elastic forces
                    call this%matrix12%matrix12_obj_stiffness(q_2(indices12a_q), q_1(indices12a_q), this%q_0(indices12a_q), q_2(indices12a_q_b), q_1(indices12a_q_b), this%q_0(indices12a_q_b), v_2(indices12a_v_a), v_1(indices12a_v_a), v_2(indices12a_v_b), v_1(indices12a_v_b), this%number_matrix12, i)
                    !< Calling subroutine to calculate the dynamic forces (inertia) 
                    call this%matrix12%matrix12_built_mass(i)
                    call this%matrix12%matrix12_mass(this%deltat, q_2(indices12a_q), q_1(indices12a_q), this%q_0(indices12a_q), q_2(indices12a_q_b), q_1(indices12a_q_b), this%q_0(indices12a_q_b), v_2(indices12a_v_a), v_1(indices12a_v_a), v_2(indices12a_v_b), v_1(indices12a_v_b), this%number_matrix12, i)
                    !< Actualization elastic forces 
                    this%fqint(indices_q_total)  = this%fqint(indices_q_total) + this%matrix12%matrix12_object(i)%fqint_obj ! fqint_obj fasst Knoten A und B zusammen
                    !< Actualization of dynamic forces
                    this%fqdyn(indices_q_total)  = this%fqdyn(indices_q_total) + this%matrix12%matrix12_object(i)%fqdyn_AB
                    !< Actualization of the forces associated with velocities 
                    this%fv(indices12v_q_total)  = this%fv(indices12v_q_total) + this%matrix12%matrix12_object(i)%fv_AB          
                    
                    !< Actualization iteration matrix
                    !< kqq (resulting from stiffnes terms) 
                    call sparse_smatrix%assemblyIndexed(indices_q_total, indices_q_total, this%matrix12%matrix12_object(i)%kqq_obj) 
                    !< kvv (resulting from inertia terms) 
                    call sparse_smatrix%assemblyIndexed(indicesv(indices12v_q_total), indicesv(indices12v_q_total), this%matrix12%matrix12_object(i)%kvv)
                    !< kqv (resulting from inertia terms)
                    call sparse_smatrix%assemblyIndexed(indices_q_total, indicesv(indices12v_q_total), this%matrix12%matrix12_object(i)%kqv)
                    !< kvq (resulting from inertia terms)
                    call sparse_smatrix%assemblyIndexed(indicesv(indices12v_q_total), indices_q_total, this%matrix12%matrix12_object(i)%kvq)              

                    this%penergy = this%penergy + this%matrix12%matrix12_object(i)%penergy  ! Hier noch aus mass erg\E4nzen
    
                end if
        end select
                
        end do
    end if 


! damping12 element
    if (allocated(this%damping12)) then
      do i = 1, this%ndamping12
        indices12a_q = this%nodes12(this%damping12(i)%c%node)%coordinates
        indices12a_v = this%nodes12(this%damping12(i)%c%node)%velocity
        
        this%damping12(i)%c%deltat = this%deltat
        
        call this%damping12(i)%c%damping12_f(v_2(indices12a_v), v_1(indices12a_v), q_2(indices12a_q), q_1(indices12a_q), this%q_0(indices12a_q))
        this%fqint(indicesq(indices12a_q)) = this%fqint(indicesq(indices12a_q)) + this%damping12(i)%c%fqint

        call this%damping12(i)%c%damping12_k(v_2(indices12a_v), v_1(indices12a_v), q_2(indices12a_q), q_1(indices12a_q), this%q_0(indices12a_q))
        call sparse_smatrix%assemblyIndexed(indicesq(indices12a_q), indicesq(indices12a_q), this%damping12(i)%c%kqq)
        call sparse_smatrix%assemblyIndexed(indicesq(indices12a_q), indicesv(indices12a_v), this%damping12(i)%c%kqv)

        this%penergy = this%penergy + this%damping12(i)%c%penergy
      end do
    end if

    return
  end subroutine modelfk

!> @brief update model_structure constraints and constraint jacobian matrix
!
!> Assembles jacobian matrix directly into the given iteration matrix
!
!> @param[in] qtn1 Current geometry configuration
!> @param[in] qtn0 last converged geometry configuration
!> @param[inout] smatrix The iteration matrix
!> @param[in] indicesg Constraint indices for assembly
!> @param[in] indicesq Positional indices for assembly
!> @param[in] lambda Lagrange multiplier vector
subroutine modelgdg(this, qtn1, qtn, sparse_smatrix, indicesg, indicesq, lambda, qopt)

    implicit none

    class(model_structure), intent(inout) :: this
    type(sparse_matrix),  intent(inout)  :: sparse_smatrix

    logical :: stiffness, boolean_q, boolean_v

    integer, allocatable, dimension(:),  intent(in) :: indicesg, indicesq
    integer :: indices12a_q(12), indices12b_q(12), indices6a(6), indices6b(6)
    integer :: constraint_indicesg(24), constraint_indicesq(24)
    integer :: i, j, rank, irank, coordinates

    real(kind = 8), intent(in) :: qtn1(:), qtn(:)
    real(kind = 8), intent(in), optional  :: qopt(:)
    real(kind = 8), optional,  intent(in) :: lambda(:)

    real(kind = 8), allocatable :: constraint_g(:), constraint_dgq(:,:), constraint_dgv(:,:), constraint_k(:,:), constraint_lambda(:)
    real(kind = 8), allocatable :: qtn05(:), qtemp(:)

    real(kind = 8), dimension(12) :: temp_qtn1_12a, temp_qtn1_12b, temp_q0_12a, temp_q0_12b, temp_qtemp_12a, temp_qtn05_12a, temp_qtn05_12b

! HACK: dummy parameter to make PARDISO work
	  real(kind = 8), allocatable :: dummy_mat(:,:)

! initialize to small value
    allocate(dummy_mat(1,1))
    dummy_mat = tiny(0.0d0)

! HACK: zeros to the main diagonal to make PARDISO work
    if (present(lambda)) then
      do i = 1, size(indicesg)
        call sparse_smatrix%assemblyIndexed(indicesg(i:i), indicesg(i:i), dummy_mat)
      end do
    endif

    allocate(constraint_g(24))
    allocate(constraint_dgq(24,24))
    allocate(constraint_dgv(24,24))
    allocate(constraint_k(24,24))
    allocate(constraint_lambda(24))
    allocate(qtn05(size(qtn1)))

    this%g(:) = 0.0d0

    allocate(qtemp(size(qtn)))
    if (present(qopt)) then
      qtemp(:) = qopt(:)
    else
      qtemp(:) = qtn(:)
    end if

    irank = 0
    if (allocated(this%constraints12) .and. (size(this%constraints12) > 0)) then

      do i = 1, size(this%constraints12)

        rank        = this%constraints12(i)%c%rankcount
        coordinates = this%constraints12(i)%c%coordinates
        stiffness   = this%constraints12(i)%c%stiffness
        boolean_q   = this%constraints12(i)%c%boolean_q
        boolean_v   = this%constraints12(i)%c%boolean_v

        this%constraints12(i)%c%deltat = this%deltat

! Indices for constraints
        constraint_indicesg = 1
        do j = 1, rank
          constraint_indicesg(j) = irank + j
        end do

! getting indices for coordinates and velocity for node a and b
        indices12a_q = this%nodes12(this%constraints12(i)%c%nodes(1))%coordinates
        indices12b_q(:) = 1

        if (coordinates > 1) then
          indices12b_q = this%nodes12(this%constraints12(i)%c%nodes(2))%coordinates
       endif

! indices for the constraint coordinate indices of node a and node b
        constraint_indicesq(1:24) = [indices12a_q, indices12b_q]

! Compute constraint equation to calculate later the constraint forces

        constraint_g(:) = 0.0d0
        constraint_dgq(:,:)  = 0.0d0
        constraint_k(:,:) = 0.0d0

        if (boolean_v .eqv. .FALSE.) then ! holonomic constraints

          ! ghol and H at tn+1
          !call this%constraints12(i)%c%gconstraint12(qtn1(indices12a_q), qtn1(indices12b_q), this%q_0(indices12a_q), this%q_0(indices12b_q), constraint_g, qtemp(indices12a_q))
          temp_qtn1_12a = qtn1(indices12a_q)
          temp_qtn1_12b = qtn1(indices12b_q)
          temp_q0_12a = this%q_0(indices12a_q)
          temp_q0_12b = this%q_0(indices12b_q)
          temp_qtemp_12a = qtemp(indices12a_q)
          call this%constraints12(i)%c%gconstraint12(temp_qtn1_12a, temp_qtn1_12b, temp_q0_12a, temp_q0_12b, constraint_g, temp_qtemp_12a)
          this%g(constraint_indicesg(1:rank)) = constraint_g(1:rank)

          call this%constraints12(i)%c%dgqconstraint12(temp_qtn1_12a, temp_qtn1_12b, temp_q0_12a, temp_q0_12b, constraint_dgq, temp_qtemp_12a)
          !call this%constraints12(i)%c%dgqconstraint12(qtn1(indices12a_q), qtn1(indices12b_q), this%q_0(indices12a_q), this%q_0(indices12b_q), constraint_dgq, qtemp(indices12a_q))
          call sparse_smatrix%assemblyIndexed(indicesg(constraint_indicesg(1:rank)), constraint_indicesq(1:coordinates*12),constraint_dgq(1:rank, 1:coordinates*12))
          
          ! H^T at tn+1/2
          constraint_dgq = 0.0d0
          qtn05 = ( qtn1 + qtn ) * 0.5d0
          !call this%constraints12(i)%c%dgqconstraint12(qtn05(indices12a_q), qtn05(indices12b_q), this%q_0(indices12a_q), this%q_0(indices12b_q), constraint_dgq, qtemp(indices12a_q))
          temp_qtn05_12a = qtn05(indices12a_q)
          temp_qtn05_12b = qtn05(indices12b_q)
          call this%constraints12(i)%c%dgqconstraint12(temp_qtn05_12a, temp_qtn05_12b, temp_q0_12a, temp_q0_12b, constraint_dgq, temp_qtemp_12a)      
          call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:coordinates*12), indicesg(constraint_indicesg(1:rank) ), transpose(constraint_dgq(1:rank, 1:coordinates*12)))

          ! Calculate d^2g/dqdq at tn+1/2
          if (stiffness) then
            constraint_lambda = lambda(constraint_indicesg(1:rank))
            call this%constraints12(i)%c%kconstraint12(constraint_lambda, constraint_k, qtemp(indices12a_q))
            call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:coordinates*12), constraint_indicesq(1:coordinates*12), 0.5d0 * constraint_k(1:coordinates*12, 1:coordinates*12) )
          end if

        elseif (boolean_v .eqv. .TRUE.) then ! nonholonomic constraints

          ! gnhol and G1, G2 at tn+1/2
          qtn05 = ( qtn1 + qtn ) * 0.5d0
          call this%constraints12(i)%c%gconstraint12(qtn05(indices12a_q), qtn05(indices12b_q), qtn(indices12a_q), qtn(indices12b_q), constraint_g, qtn1(indices12a_q))
          this%g(constraint_indicesg(1:rank)) = constraint_g(1:rank)

          call this%constraints12(i)%c%dgvconstraint12(qtn05(indices12a_q), qtn05(indices12b_q), qtn(indices12a_q), qtn(indices12b_q), constraint_dgv, qtn1(indices12a_q))
          call sparse_smatrix%assemblyIndexed(indicesg(constraint_indicesg(1:rank) ), constraint_indicesq(1:coordinates*12), constraint_dgv(1:rank, 1:coordinates*12))

          ! G^T (=G2^T) at tn+1/2
          constraint_dgq = 0.0d0
          call this%constraints12(i)%c%dgqconstraint12(qtn05(indices12a_q), qtn05(indices12b_q), qtn(indices12a_q), qtn(indices12b_q), constraint_dgq, qtn1(indices12a_q))
          call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:coordinates*12), indicesg(constraint_indicesg(1:rank) ), transpose(constraint_dgq(1:rank, 1:coordinates*12)))

          ! Calculate d^2g/dqdq at tn+1/2
          if (stiffness) then
            constraint_lambda = lambda(constraint_indicesg(1:rank))
            call this%constraints12(i)%c%kconstraint12(constraint_lambda, constraint_k, qtemp(indices12a_q))
            call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:coordinates*12), constraint_indicesq(1:coordinates*12), 0.5d0 * constraint_k(1:coordinates*12, 1:coordinates*12) )
          end if

        end if

        irank = irank+rank
      end do
    end if

! node 6
    if (allocated(this%constraints6) .and. (size(this%constraints6) > 0)) then

      do i = 1, size(this%constraints6)

        rank = this%constraints6(i)%c%rankcount
        coordinates = this%constraints6(i)%c%coordinates
        stiffness = this%constraints6(i)%c%stiffness

        constraint_indicesg = 1
        do j = 1, rank
          constraint_indicesg(j) = irank+j
        end do
        indices6a = this%nodes6(this%constraints6(i)%c%nodes(1))%coordinates

        if (this%constraints6(i)%c%nodes(2) > 0) then
          indices6b = this%nodes6(this%constraints6(i)%c%nodes(2))%coordinates
        else
          indices6b = indices6a
        endif
        constraint_indicesq(1:12) = [indices6a, indices6b]

        ! g at tn+1
        constraint_g = 0.0d0
        call this%constraints6(i)%c%gconstraint6(qtn1(indices6a), qtn1(indices6b), this%q_0(indices6a), this%q_0(indices6b), constraint_g, qtemp(indices6a) )
        this%g(constraint_indicesg(1:rank)) = constraint_g(1:rank)

        ! H at tn+1
        constraint_dgq = 0.0d0
        call this%constraints6(i)%c%dgconstraint6(qtn1(indices6a), qtn1(indices6b), this%q_0(indices6a), this%q_0(indices6b), constraint_dgq, qtemp(indices6a))
        call sparse_smatrix%assemblyIndexed(indicesg(constraint_indicesg(1:rank) ), constraint_indicesq(1:coordinates*6), constraint_dgq(1:rank, 1:coordinates*6))

        ! H^T at tn+1/2
        qtn05 = ( qtn1 + qtn ) * 0.5d0
        constraint_dgq = 0.0d0
        call this%constraints6(i)%c%dgconstraint6(qtn05(indices6a), qtn05(indices6b), this%q_0(indices6a), this%q_0(indices6b), constraint_dgq, qtemp(indices6a))
        call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:coordinates*6), indicesg(constraint_indicesg(1:rank) ), transpose(constraint_dgq(1:rank, 1:coordinates*6)))

        if (stiffness) then
          constraint_k = 0.0d0
          constraint_lambda = lambda(constraint_indicesg(1:rank))
          call this%constraints6(i)%c%kconstraint6(constraint_lambda, constraint_k)
          call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:coordinates*6), constraint_indicesq(1:coordinates*6), 0.5d0 * constraint_k(1:coordinates*6, 1:coordinates*6) )
        end if

        irank = irank+rank

       end do

    end if

! node 6 to 12
    if (allocated(this%constraints6to12) .and. (size(this%constraints6to12) > 0)) then

      do i = 1, size(this%constraints6to12)

        rank =      this%constraints6to12(i)%c%rankcount
        stiffness = this%constraints6to12(i)%c%stiffness

        do j = 1, rank
          constraint_indicesg(j) = irank+j
        end do
        indices12a_q = this%nodes12(this%constraints6to12(i)%c%nodes(1))%coordinates
        indices6b  = this%nodes6 (this%constraints6to12(i)%c%nodes(2))%coordinates
        constraint_indicesq(1:18) = [indices12a_q, indices6b]

        ! g at tn+1
        constraint_g = 0.0d0
        call this%constraints6to12(i)%c%gconstraint6to12(qtn1(indices12a_q), qtn1(indices6b), this%q_0(indices12a_q), this%q_0(indices6b), constraint_g )
        this%g (constraint_indicesg(1:rank)) = constraint_g(1:rank)

        ! H at tn+1
        constraint_dgq = 0.0d0
        call this%constraints6to12(i)%c%dgconstraint6to12(qtn1(indices12a_q), qtn1(indices6b), this%q_0(indices12a_q), this%q_0(indices6b), constraint_dgq)
        call sparse_smatrix%assemblyIndexed(indicesg(constraint_indicesg(1:rank) ), constraint_indicesq(1:18), constraint_dgq(1:rank, 1:18))

        ! H^T at tn+1/2
        qtn05 = ( qtn1 + qtn ) * 0.5d0
        constraint_dgq = 0.0d0
        call this%constraints6to12(i)%c%dgconstraint6to12(qtn05(indices12a_q), qtn05(indices6b), this%q_0(indices12a_q), this%q_0(indices6b), constraint_dgq)
        call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:18), indicesg(constraint_indicesg(1:rank) ), transpose(constraint_dgq(1:rank, 1:18)))

        if (stiffness) then
          constraint_k = 0.0d0
          constraint_lambda = lambda(constraint_indicesg(1:rank))
          call this%constraints6to12(i)%c%kconstraint6to12(constraint_lambda, constraint_k)
          call sparse_smatrix%assemblyIndexed(constraint_indicesq(1:18), constraint_indicesq(1:18), 0.5d0 * constraint_k(1:18, 1:18) )
        end if

        irank = irank+rank

      end do

    end if

    return

end subroutine modelgdg

!< Model null space
subroutine modelndgq(this, qtn1, qtn, sparse_nullspace, qopt)

    use my_math_structure, only: eye
    use sparse_function

    implicit none

    class(model_structure), intent(inout) :: this
    class(single_sparse_matrix), intent(inout) :: sparse_nullspace
    integer :: indices12a_q(12), indices12b_q(12), indices6a_q(6), indices6b_q(6)
    integer :: i, j, k, ie, r1
    integer :: c12, c6
    integer :: n_nodes12, n_nodes6
    integer :: rank, node1, node2, i_c, rank_j, constraint_j, coordinates_j, n_q
    integer, allocatable :: indices_c(:), indices_q_n(:), indices_temp_q(:), indices_dgq_local(:)
    integer :: indices_q_n_add(6)

    real(kind = 8), intent(in) :: qtn1(:), qtn(:)
    real(kind = 8), intent(in), optional :: qopt(:)

    real(kind = 8), allocatable :: constraint_dgq(:,:), constraint_ndgq(:,:)
    real(kind = 8), allocatable :: H_i(:,:)
    real(kind = 8), allocatable :: qtemp(:)

    allocate(constraint_dgq(24,24))

    allocate(qtemp(size(qtn1)))
    if (present(qopt)) then
      qtemp(:) = qopt(:)
    else
      qtemp(:) = qtn(:)
    end if
    
    if (allocated(this%null_space%list)) then
      ! Loop over the costraint nullspace list
      do i = 1, size(this%null_space%list)

            n_nodes12 = this%null_space%list(i)%n_nodes12
            n_nodes6  = this%null_space%list(i)%n_nodes6

            rank = this%null_space%list(i)%rank
            n_q  =  size(this%null_space%list(i)%nodes12)*12 + size(this%null_space%list(i)%nodes6)*6

            if (rank /= n_q) then
                  ! initialize constraint matrix
                  allocate(H_i(rank,n_q))
                  H_i = 0.0d0
                  i_c = 0
                  do j = 1,size(this%null_space%list(i)%constraints)

                        constraint_j = this%null_space%constraint_list(this%null_space%list(i)%constraints(j),5)

                        c12 = count(this%null_space%constraint_list(this%null_space%list(i)%constraints(j),1:2) .GT. 0)
                        c6  = count(this%null_space%constraint_list(this%null_space%list(i)%constraints(j),3:4) .GT. 0)

                        if (c12 > 0 .AND. c6 > 0) then
                        ! constraint6to12 - for beams, rigid bodies and shell nodes
                              coordinates_j = this%constraints6to12(constraint_j)%c%coordinates
                              rank_j        = this%constraints6to12(constraint_j)%c%rankcount

                              node1 = this%constraints6to12(constraint_j)%c%nodes(1)
                              node2 = this%constraints6to12(constraint_j)%c%nodes(2)

                              ! global node indices
                              indices12a_q = this%nodes12(node1)%coordinates
                              indices6b_q  = this%nodes6(node2)%coordinates

                              ! Calculate dg/dq
                              constraint_dgq = 0.0d0
                              call this%constraints6to12(constraint_j)%c%dgconstraint6to12(qtn1(indices12a_q), qtn1(indices6b_q), this%q_0(indices12a_q), this%q_0(indices6b_q), constraint_dgq)

                              ! local indices for this subroutine
                              allocate(indices_temp_q(18))
                              indices_temp_q(:) = 1
                              do k = 1,n_nodes12 !size(this%null_space%list(i)%nodes12)
                                if (node1 == this%null_space%list(i)%nodes12(k)) then
                                  indices_temp_q(1:12) = this%null_space%list(i)%node(k)%indices_temp_q(:)
                                  exit
                                end if
                              end do
                              do k = 1,n_nodes6 !size(this%null_space%list(i)%nodes6)
                                if (node2 == this%null_space%list(i)%nodes6(k)) then
                                  indices_temp_q(13:18) = this%null_space%list(i)%node(size(this%null_space%list(i)%nodes12)+k)%indices_temp_q(:)
                                  exit
                                end if
                              end do

                              ! local indices for constraint matrix dgq
                              allocate(indices_dgq_local(18))
                              indices_dgq_local(1:18) = (/ (r1, r1 = 1,18) /)

                        else if (c12 > 0 .AND. c6 == 0) then
                              ! constraint12 - for beams and rigid body nodes
                              node1 = this%constraints12(constraint_j)%c%nodes(1)
                              node2 = this%constraints12(constraint_j)%c%nodes(2)

                              coordinates_j = this%constraints12(constraint_j)%c%coordinates
                              rank_j        = this%constraints12(constraint_j)%c%rankcount

                              ! global node indices
                              indices12a_q = this%nodes12(node1)%coordinates
                              indices12b_q = 1
                              if (coordinates_j /= 1) then
                                indices12b_q(:) = this%nodes12(node2)%coordinates
                              end if

                              constraint_dgq = 0.0d0
                              call this%constraints12(constraint_j)%c%dgqconstraint12(qtn1(indices12a_q), qtn1(indices12b_q), this%q_0(indices12a_q), this%q_0(indices12b_q), constraint_dgq, qtemp(indices12a_q))

                              ! local indices for this subroutine
                              allocate(indices_temp_q(coordinates_j*12))
                              indices_temp_q(:) = 1
                            
                              do k = 1,n_nodes12 !size(this%null_space%list(i)%nodes12)
                                if (node1 == this%null_space%list(i)%nodes12(k)) then
                                  indices_temp_q(1:12) = this%null_space%list(i)%node(k)%indices_temp_q(:)
                                end if
                              end do
                            
                              if (coordinates_j == 2) then
                                do k = 1,n_nodes12 ! size(this%null_space%list(i)%nodes12)
                                  if (node2 == this%null_space%list(i)%nodes12(k)) then
                                    indices_temp_q(13:24) = this%null_space%list(i)%node(k)%indices_temp_q(:)
                                  end if
                                end do
                              end if
                              ! local indices for constraint matrix dgq
                              allocate(indices_dgq_local(12*coordinates_j))
			                        indices_dgq_local(1:12*coordinates_j) = (/ (r1, r1 = 1,12*coordinates_j) /)

                        elseif (c12 == 0 .AND. c6 > 0) then
                              ! constraints6 - for shelle nodes
                              node1 = this%constraints6(constraint_j)%c%nodes(1)
                              node2 = this%constraints6(constraint_j)%c%nodes(2)

                              coordinates_j = this%constraints6(constraint_j)%c%coordinates
                              rank_j        = this%constraints6(constraint_j)%c%rankcount

                              ! global node indices
                              indices6a_q = this%nodes6(node1)%coordinates
                              indices6b_q = 1
                              if (coordinates_j /= 1) then
                                indices6b_q(:) = this%nodes6(node2)%coordinates
                              else
                                indices6b_q = indices6a_q
                              end if

                              constraint_dgq = 0.0d0
                              call this%constraints6(constraint_j)%c%dgconstraint6(qtn1(indices6a_q), qtn1(indices6b_q), this%q_0(indices6a_q), this%q_0(indices6b_q), constraint_dgq, qtemp(indices6a_q))

                              ! local indices for this subroutine
                              allocate(indices_temp_q(coordinates_j*6))
                              indices_temp_q(:) = 1
                              do k = 1,n_nodes6 !size(this%null_space%list(i)%nodes6)
                                if (node1 == this%null_space%list(i)%nodes6(k)) then
                                  indices_temp_q(1:6) = this%null_space%list(i)%node(size(this%null_space%list(i)%nodes12)+k)%indices_temp_q(:)
                                end if
                              end do
                              if (coordinates_j == 2) then
                                do k = 1,n_nodes6 !size(this%null_space%list(i)%nodes6)
                                  if (node2 == this%null_space%list(i)%nodes6(k)) then
                                    indices_temp_q(7:12) = this%null_space%list(i)%node(size(this%null_space%list(i)%nodes12)+k)%indices_temp_q(:)
                                  end if
                                end do
                              end if
                              ! local indices for constraint matrix dgq
                              allocate(indices_dgq_local(6*coordinates_j))
			                        indices_dgq_local(1:6*coordinates_j) = (/ (r1, r1 = 1,6*coordinates_j) /)
                        end if

                        ! creating constraint indices, i_c is sequential +1
                        allocate(indices_c(rank_j))
                        indices_c = 0
                        do k = 1, rank_j
                          i_c = i_c + 1
                          indices_c(k) = i_c
                        end do

                        ! calculation of constraint matrix
                        H_i(indices_c, indices_temp_q) = constraint_dgq(1:rank_j,indices_dgq_local)

                        deallocate(indices_c)
                        deallocate(indices_temp_q)
                        deallocate(indices_dgq_local)
                  end do

                  ! getting global indices of the constraint nodes
                  allocate(indices_q_n(n_q))
                  indices_q_n = 0
                  ie = 0
                  do k = 1,n_nodes12 !size(this%null_space%list(i)%nodes12)
                    indices_q_n(ie+1:ie+12) = this%nodes12(this%null_space%list(i)%nodes12(k))%coordinates
                    ie = ie + 12
                  end do

                  do k = 1,n_nodes6 !size(this%null_space%list(i)%nodes6)
                    indices_q_n(ie+1:ie+6) = this%nodes6(this%null_space%list(i)%nodes6(k))%coordinates_q
                    ie = ie + 6
                  end do

                  ! calculate null space using lapack function: get_nullspace_mod(N_matrix,matrix,nq,nc)
                  allocate(constraint_ndgq(n_q,n_q-rank))
                  constraint_ndgq = 0.0d0
                  call get_nullspace_mod(constraint_ndgq, H_i, n_q, rank)

                  ! into sparse null space matrix:  assemblyIndexed(row_indices, col_indices, m)
                  call create_rowcol_format(sparse_nullspace, indices_q_n, this%null_space%list(i)%indicesn,constraint_ndgq)

                  ! Deallocating variables
                  deallocate(H_i)
                  deallocate(indices_q_n)
                  deallocate(constraint_ndgq)

            end if
      end do
    end if

    ! This is only for node6, because for node12 all nodes are always constrained due to internal constraint
    if (allocated(this%null_space%addlist6)) then
      do i = 1,size(this%null_space%addlist6)
        ! into sparse null space matrix:  assemblyIndexed(row_indices, col_indices, m)
        indices_q_n_add = this%nodes6(this%null_space%addlist6(i))%coordinates_q
        call create_rowcol_format(sparse_nullspace, indices_q_n_add, this%null_space%ns_add6(i)%indicesn,eye(6))
      end do
    end if

    return
  end subroutine modelndgq

! this function compute the total rank of the constraints.
subroutine constraintsrankcounter(this)
    class(model_structure), intent(inout) :: this
    integer :: i, irank

    irank = 0
    if (allocated(this%constraints6) .and. (size(this%constraints6) > 0)) then
       do i = 1, size(this%constraints6)
          irank = irank+this%constraints6(i)%c%rankcount
       end do
    end if
    this%trconstraints6 = irank
    
    irank = 0
    if (allocated(this%constraints12) .and. (size(this%constraints12) > 0)) then
       do i = 1, size(this%constraints12)
          irank = irank+this%constraints12(i)%c%rankcount
       end do
    end if
    this%trconstraints12 = irank
    
    irank = 0
    if (allocated(this%constraints6to12) .and. (size(this%constraints6to12) > 0)) then
       do i = 1, size(this%constraints6to12)
          irank = irank+this%constraints6to12(i)%c%rankcount
       end do
    end if
    this%trconstraints6to12 = irank
    
    this%trconstraints = this%trconstraints6+this%trconstraints12+this%trconstraints6to12

    return
  end subroutine constraintsrankcounter

subroutine modelloads(this, time, ml, sl)
    implicit none

    class(model_structure), intent(in) :: this
    real(kind = 8) :: time
    real(kind = 8), allocatable, dimension(:), intent(inout) :: ml, sl
    integer :: i, i6

    ml(:) = 0.0d0
    sl(:) = 0.0d0

    if(allocated(this%loads12) .and. this%nloads12 > 0) then
       do i = 1, this%nloads12
          i6 = 6*(this%loads12(i)%node-1)
          ml(i6+1:i6+6) = ml(i6+1:i6+6) + this%loads12(i)%material*this%loadamplitudes12(i)%amplitude(time)
          sl(i6+1:i6+6) = sl(i6+1:i6+6) + this%loads12(i)%spatial *this%loadamplitudes12(i)%amplitude(time)
       end do
    end if

    if(allocated(this%loads6) .and. this%nloads6 > 0) then
       do i = 1, this%nloads6
          i6 = 6*(this%loads6(i)%node-1)+6*this%nnodes12
          ml(i6+3     ) = ml(i6+3     ) + this%loads6(i)%material*this%loadamplitudes6(i)%amplitude(time)
          sl(i6+1:i6+3) = sl(i6+1:i6+3) + this%loads6(i)%spatial *this%loadamplitudes6(i)%amplitude(time)
       end do
    end if

    return
  end subroutine modelloads

subroutine modelboundaries(this, time)
    implicit none

    class(model_structure), intent(inout) :: this
    real(kind = 8), intent(in) :: time
    real(kind = 8) :: amplitude, amplitude0
    integer :: i

    if(allocated(this%constraints12)) then

       do i = 1, size(this%constraints12)
         this%constraints12(i)%c%amplitude_translation     = 0.0d0
         this%constraints12(i)%c%amplitude_translation_vel = 0.0d0
         this%constraints12(i)%c%amplitude_rotation        = 0.0d0
         this%constraints12(i)%c%amplitude_rotation_vel    = 0.0d0
       end do

    end if

    if(allocated(this%boundaries12)) then

       do i = 1, size(this%boundaries12)
         ! calculate amplitude of manipulation
         amplitude = this%boundaryamplitudes12(i)%amplitude(time)
         if (time == 0.0d0) then
            amplitude0 = 0.0d0
         else
            amplitude0 = this%boundaryamplitudes12(i)%amplitude(time-this%deltat)
         end if

         this%constraints12(this%boundaries12(i)%constraint12)%c%amplitude_translation     = amplitude
         this%constraints12(this%boundaries12(i)%constraint12)%c%amplitude_translation_vel = amplitude
         this%constraints12(this%boundaries12(i)%constraint12)%c%amplitude_rotation        = (amplitude-amplitude0)
         this%constraints12(this%boundaries12(i)%constraint12)%c%amplitude_rotation_vel    = amplitude

       end do

    end if

    ! boundary manipulation for 6 DOF nodes
    if(allocated(this%constraints6)) then

       do i = 1, size(this%constraints6)
         this%constraints6(i)%c%amplitude_translation     = 0.0d0
         this%constraints6(i)%c%amplitude_translation_vel = 0.0d0
         this%constraints6(i)%c%amplitude_rotation        = 0.0d0
         this%constraints6(i)%c%amplitude_rotation_vel    = 0.0d0
       end do

    end if

    if(allocated(this%boundaries6)) then

       do i = 1, size(this%boundaries6)
         ! calculate amplitude of manipulation
         amplitude = this%boundaryamplitudes6(i)%amplitude(time)
         if (time == 0.0d0) then
            amplitude0 = 0.0d0
         else
            amplitude0 = this%boundaryamplitudes6(i)%amplitude(time-this%deltat)
         end if

         this%constraints6(this%boundaries6(i)%constraint6)%c%amplitude_translation     = amplitude
         this%constraints6(this%boundaries6(i)%constraint6)%c%amplitude_translation_vel = amplitude
         this%constraints6(this%boundaries6(i)%constraint6)%c%amplitude_rotation        = amplitude
         this%constraints6(this%boundaries6(i)%constraint6)%c%amplitude_rotation_vel    = amplitude
                  
       end do

    end if

    return

  end subroutine modelboundaries

!< Subroutine for reading input files
  subroutine readinginputs(this)

    use my_math_structure, only: vec6mat9symm, vec10mat16symm
    use my_materials, only: orthotropicmaterial

    implicit none

    class(model_structure), intent(inout) :: this

    integer :: io_error
    integer :: i, j, k, l
    integer, parameter :: linejump = 3

    integer :: nbodyproperties, nbeamproperties, nshellmaterials
    real(kind = 8), allocatable :: matpropbeam(:,:)
    real(kind = 8), allocatable :: cmass10v(:, :), cmatv(:, :), density(:), alphashelldiss(:,:), alphabeamdiss(:,:)
    real(kind = 8) :: celast(6,6)
    integer :: nshellsequences
    integer, allocatable :: shellsequences(:, :)
    real(kind = 8), allocatable :: anglesequences(:, :), thicknesssequences(:, :)

    character(50) :: constraintsort, loadsort, stepsort, filename_boundary12, filename_boundary6, spring12_sort, damping12_sort
    character(100) :: strfilename, stroldfilename, strtype, strtemp
    
    integer :: nStepCount
    type(constraint) :: temp_constraint

    real(kind = 8) ::  spring12dir(3), damping12dir(3)
    integer :: spring12node, spring12property, damping12node, damping12property
    type(condensed_12_curve), allocatable :: spring12properties(:), damping12properties(:)

    integer :: npostproc
    integer :: nnodes
    integer, allocatable :: postprocid(:)
    integer :: flag_property_type, io
    
!< reading rigid bodies
    this%nbodies = 0
    nbodyproperties = 0
    open(unit = 1000, file = 'rigidbodyinput.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then ! if the file rigidbodyinput exists, then read the information.
       do i = 1, linejump
          read(1000, *)
       end do
       read(1000, *) this%nbodies, nbodyproperties
       allocate(this%bodies(this%nbodies))
       allocate(cmass10v(nbodyproperties, 10))
       do i = 1, this%nbodies
          do j = 1, linejump
             read(1000, *)
          end do
          read(1000, *) this%bodies(i)%node%q_0(:)
          do j = 1, linejump
             read(1000, *)
          end do
          read(1000, *) this%bodies(i)%property
       end do
       do i = 1, nbodyproperties
          do j = 1, linejump
             read(1000, *)
          end do
          read(1000, *) cmass10v(i, :)
       end do
       ! accomodating information
       do i = 1, this%nbodies
          j = this%bodies(i)%property
          this%bodies(i)%cmass = vec10mat16symm(cmass10v(j, :))
       end do
       deallocate(cmass10v)
    end if
    close(unit = 1000)
    
!< reading beams
    nnodes = this%nbodies
    this%nbeams = 0
    nbeamproperties = 0
    open(unit = 2000, file = 'beaminput.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file beaminput exists, then read the information.
       do i = 1, linejump
          read(2000, *)
       end do

       ! Number of beams (body) and number of material/section properties
       flag_property_type = 0
       read(2000, *) this%nbeams, nbeamproperties, strtemp
       
       ! set input property type
       if (trim(adjustl(strtemp)) .eq. '!!') then
         backspace(2000)
       else
         read(strtemp, '(i1)') flag_property_type
       end if
       
       allocate(this%beams(this%nbeams))

       do i = 1, this%nbeams
          do j = 1, linejump
             read(2000, *)
      end do

          ! number of beam elements of body beam and number of nodes of beam elements
          read(2000, *) this%beams(i)%nnodes, this%beams(i)%nelements

          allocate(this%beams(i)%nodes   (this%beams(i)%nnodes   ))
          allocate(this%beams(i)%elements(this%beams(i)%nelements))

          do j = 1, linejump
             read(2000, *)
          end do

          ! reading nodes
          do j = 1, this%beams(i)%nnodes
             read(2000, *) this%beams(i)%nodes(j)%q_0(:)
          end do

          do j = 1, linejump
             read(2000, *)
          end do

          ! reading elements
          do j = 1, this%beams(i)%nelements
             read(2000, *) this%beams(i)%elements(j)%connectivity, this%beams(i)%elements(j)%property
          end do
          
          nnodes = nnodes + this%beams(i)%nnodes
       end do
       
       if (flag_property_type .eq. 0) then
          allocate(matpropbeam(nbeamproperties,17))
       elseif (flag_property_type .eq. 1) then
          allocate(matpropbeam(nbeamproperties,21+6))
       end if
       allocate(alphabeamdiss(nbeamproperties,2))       
       matpropbeam   = 0.0d0
       alphabeamdiss = 0.0d0
       
       ! reading beam properties
       do i = 1, nbeamproperties
          do j = 1, linejump
             read(2000, *)
          end do

          if (flag_property_type .eq. 0) then
            read(2000, *) matpropbeam(i, 1), matpropbeam(i, 2), matpropbeam(i, 3), matpropbeam(i, 4), matpropbeam(i, 5), matpropbeam(i, 6), matpropbeam(i, 7), matpropbeam(i, 8), matpropbeam(i, 9), matpropbeam(i, 10), matpropbeam(i, 11)   ! EA, GA1, GA2, EI1, EI2, GIP, ES1, ES2, GS1, GS2, EI12
            read(2000, *) matpropbeam(i, 12), matpropbeam(i, 13), matpropbeam(i, 14), matpropbeam(i, 15), matpropbeam(i, 16), matpropbeam(i, 17)  ! rhoA3, rhoI1, rhoI2, rhoS1, rhoS2, rhoI12
          elseif (flag_property_type .eq. 1) then
            read(2000, *) matpropbeam(i, 1), matpropbeam(i, 2), matpropbeam(i, 3), matpropbeam(i, 4), matpropbeam(i, 5), matpropbeam(i, 6), matpropbeam(i, 7), matpropbeam(i, 8), matpropbeam(i, 9), matpropbeam(i, 10), matpropbeam(i, 11), matpropbeam(i, 12), matpropbeam(i, 13), matpropbeam(i, 14), matpropbeam(i, 15), matpropbeam(i, 16), matpropbeam(i, 17), matpropbeam(i, 18), matpropbeam(i, 19), matpropbeam(i, 20), matpropbeam(i, 21)   ! all entries in elasticity matrix
            read(2000, *) matpropbeam(i, 22), matpropbeam(i, 23), matpropbeam(i, 24), matpropbeam(i, 25), matpropbeam(i, 26), matpropbeam(i, 27)  ! all entries in mass matrix
          end if

          ! if (flag_property_type .eq. 0) then
          !   read(2000, *) matpropbeam(i,1:11)   ! EA, GA1, GA2, EI1, EI2, GIP, ES1, ES2, GS1, GS2, EI12
          !   read(2000, *) matpropbeam(i,12:17)  ! rhoA3, rhoI1, rhoI2, rhoS1, rhoS2, rhoI12
          ! elseif (flag_property_type .eq. 1) then
          !   read(2000, *) matpropbeam(i,1:21)   ! all entries in elasticity matrix
          !   read(2000, *) matpropbeam(i,22:27)  ! all entries in mass matrix
          ! end if
          ! diss sigma, diss velocity
          !read(2000, *) alphabeamdiss(i,:)
          read(2000, *) alphabeamdiss(i,1), alphabeamdiss(i,2)
       end do

!! accomodating information
       do i = 1, this%nbeams
          do j = 1, this%beams(i)%nelements
             k = this%beams(i)%elements(j)%property
             if (flag_property_type .eq. 0) then
               this%beams(i)%elements(j)%cgg      = 0.0d0
               this%beams(i)%elements(j)%cgg(1,1) = matpropbeam(k,2)  ! G A1
               this%beams(i)%elements(j)%cgg(2,2) = matpropbeam(k,3)  ! G A2
               this%beams(i)%elements(j)%cgg(3,3) = matpropbeam(k,1)  ! E A3

               this%beams(i)%elements(j)%ckk      =  0.0d0
               this%beams(i)%elements(j)%ckk(1,1) =  matpropbeam(k,4)  ! E I1
               this%beams(i)%elements(j)%ckk(2,2) =  matpropbeam(k,5)  ! E I2
               this%beams(i)%elements(j)%ckk(3,3) =  matpropbeam(k,6)  ! G IT
               this%beams(i)%elements(j)%ckk(1,2) = -matpropbeam(k,11) ! E I12
               this%beams(i)%elements(j)%ckk(2,1) = -matpropbeam(k,11) ! E I12

               this%beams(i)%elements(j)%cgk      =  0.0d0
               this%beams(i)%elements(j)%cgk(1,3) = -matpropbeam(k,9)  ! G S1
               this%beams(i)%elements(j)%cgk(2,3) =  matpropbeam(k,10) ! G S2
               this%beams(i)%elements(j)%cgk(3,1) =  matpropbeam(k,7)  ! E S1
               this%beams(i)%elements(j)%cgk(3,2) = -matpropbeam(k,8)  ! E S2

               this%beams(i)%elements(j)%cmass(:,:) = 0.0d0
               this%beams(i)%elements(j)%cmass(1,1) = matpropbeam(k,12)  ! density A3
               this%beams(i)%elements(j)%cmass(1,2) = matpropbeam(k,16)  ! density S2
               this%beams(i)%elements(j)%cmass(1,3) = matpropbeam(k,15)  ! density S1

               this%beams(i)%elements(j)%cmass(2,1) = matpropbeam(k,16)  ! density S2
               this%beams(i)%elements(j)%cmass(2,2) = matpropbeam(k,14)  ! density I2
               this%beams(i)%elements(j)%cmass(2,3) = matpropbeam(k,17)  ! density I12

               this%beams(i)%elements(j)%cmass(3,1) = matpropbeam(k,15)  ! density S1
               this%beams(i)%elements(j)%cmass(3,2) = matpropbeam(k,17)  ! density I12
               this%beams(i)%elements(j)%cmass(3,3) = matpropbeam(k,13)  ! density I1
             
             elseif (flag_property_type .eq. 1) then
               celast = vec21mat36symm(matpropbeam(k,1:21))
               this%beams(i)%elements(j)%cgg   = celast(1:3,1:3)
               this%beams(i)%elements(j)%ckk   = celast(4:6,4:6)
               this%beams(i)%elements(j)%cgk   = celast(1:3,4:6)
               this%beams(i)%elements(j)%cmass = vec6mat9symm(matpropbeam(k,22:27))
             end if
             
             this%beams(i)%elements(j)%alpha_s = alphabeamdiss(k,1)    ! damping (dissipation) factor regarding to stress (strain energy)
             this%beams(i)%elements(j)%alpha_v = alphabeamdiss(k,2)    ! damping (dissipation) factor regarding to velocity (kinetic energy)
          end do
       end do

       deallocate(matpropbeam)
       deallocate(alphabeamdiss)

    end if

    close(unit = 2000)

!< reading shells
    this%nshells = 0
    nshellmaterials = 0
    nshellsequences = 0
    print *, "Parsing shell element input..."
    open(unit = 3000, file = 'shellinput.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then ! if the file shellinput exists, then read the information.
       do i = 1, linejump
          read(3000, *)
       end do
       read(3000, *) this%nshells, nshellmaterials, nshellsequences
       allocate(this%shells       (this%nshells    ))
       allocate(cmatv             (nshellmaterials, 9))
       allocate(density           (nshellmaterials))
       allocate(alphashelldiss    (nshellmaterials,2))
       allocate(shellsequences    (nshellsequences, 11))
       allocate(anglesequences    (nshellsequences, 10))
       allocate(thicknesssequences(nshellsequences, 10))

       do i = 1, this%nshells
          do j = 1, linejump
             read(3000, *)
          end do
          read(3000, *) this%shells(i)%nnodes, this%shells(i)%nelements
          allocate(this%shells(i)%nodes   (this%shells(i)%nnodes   ))
          allocate(this%shells(i)%elements(this%shells(i)%nelements))
          do j = 1, linejump
             read(3000, *)
          end do
          do j = 1, this%shells(i)%nnodes
             read(3000, *) this%shells(i)%nodes(j)%q_0(:)
          end do
          do j = 1, linejump
             read(3000, *)
          end do
          do j = 1, this%shells(i)%nelements
             read(3000, *) this%shells(i)%elements(j)%connectivity, this%shells(i)%elements(j)%laminate%plysequence
          end do
          nnodes = nnodes + this%shells(i)%nnodes
       end do
       do i = 1, nshellmaterials
          do j = 1, linejump
             read(3000, *)
          end do
          read(3000, *) cmatv         (i, :)
          read(3000, *) density       (i)
          read(3000, *) alphashelldiss(i,:)
       end do
       do i = 1, nshellsequences
          do j = 1, linejump
             read(3000, *)
          end do
          read(3000, *) shellsequences    (i, 11)
          read(3000, *) shellsequences    (i, 1:shellsequences(i, 11))
          read(3000, *) thicknesssequences(i, 1:shellsequences(i, 11))
          read(3000, *) anglesequences    (i, 1:shellsequences(i, 11))
       end do
       ! accomodating information.
       do i = 1, this%nshells
          do j = 1, this%shells(i)%nelements
             k = this%shells(i)%elements(j)%laminate%plysequence
             this%shells(i)%elements(j)%laminate%nplys = shellsequences(k, 11)
             allocate(this%shells(i)%elements(j)%laminate%plys(this%shells(i)%elements(j)%laminate%nplys))
             do l = 1, this%shells(i)%elements(j)%laminate%nplys
                this%shells(i)%elements(j)%laminate%plys(l)%material  = shellsequences(k, l)
                call this%shells(i)%elements(j)%laminate%plys(l)%initialization(cmatv             (this%shells(i)%elements(j)%laminate%plys(l)%material, :), &
                                                                                density           (this%shells(i)%elements(j)%laminate%plys(l)%material), &
                                                                                thicknesssequences(k, l), &
                                                                                anglesequences    (k, l), &
                                                                                alphashelldiss    (this%shells(i)%elements(j)%laminate%plys(l)%material,1), &   ! alpha_s
                                                                                alphashelldiss    (this%shells(i)%elements(j)%laminate%plys(l)%material,2))     ! alpha_v
             end do
             allocate(this%shells(i)%elements(j)%laminate%mastermasterjacobean(this%shells(i)%elements(j)%laminate%nplys))
             allocate(this%shells(i)%elements(j)%laminate%boundaries(this%shells(i)%elements(j)%laminate%nplys+1))
             allocate(this%shells(i)%elements(j)%laminate%masterboundaries(this%shells(i)%elements(j)%laminate%nplys+1))
             call this%shells(i)%elements(j)%laminate%initialization()
          end do
       end do

       deallocate(cmatv             )
       deallocate(density           )
       deallocate(alphashelldiss    )
       deallocate(shellsequences    )
       deallocate(anglesequences    )
       deallocate(thicknesssequences)

    end if
    close(unit = 3000)

!< reading spring_12
    this%nspring12 = 0
    this%nspring12properties = 0
    print *, "Parsing springs for 12 DOF ..."
    open(unit = 3100, file = 'spring12input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file exists, then read the information.

      do i = 1, linejump
        read(3100, *)
      end do

      read(3100, *) this%nspring12, this%nspring12properties
      do i = 1, linejump
        read(3100, *)
      end do

      allocate(this%spring12(this%nspring12))
      do i = 1, this%nspring12
        read(3100, *) spring12_sort, spring12node, spring12property, spring12dir(:)

        select case (trim(adjustl(spring12_sort)))
          case ('translatory_inplane')
            allocate(this%spring12(i)%c, source = spring12_translatory_inplane() )

          case ('translatory_global')
            allocate(this%spring12(i)%c, source = spring12_translatory_global() )

          case ('translatory_local')
            allocate(this%spring12(i)%c, source = spring12_translatory_local() )

          case ('translatory_corotational')
            allocate(this%spring12(i)%c, source = spring12_translatory_corotational() )

          case default
              print *, "ERROR: didn't recognize spring12: ", spring12_sort
        end select

        allocate(this%spring12(i)%c%sort, source = trim(adjustl(spring12_sort)))
        this%spring12(i)%c%node                  = spring12node
        this%spring12(i)%c%propertynumber        = spring12property
        this%spring12(i)%c%dir(:)                = spring12dir(:)

      end do

      allocate(spring12properties(this%nspring12properties))
      do i = 1,this%nspring12properties
        do j = 1, linejump
          read(3100, *)
        end do

        read(3100, *) spring12properties(i)%npoints
        allocate(spring12properties(i)%data_points(spring12properties(i)%npoints,2))
        do j = 1,spring12properties(i)%npoints
          read(3100, *) spring12properties(i)%data_points(j,:)
        end do
      end do

      do i = 1,this%nspring12
        allocate(this%spring12(i)%c%data_points(spring12properties(this%spring12(i)%c%propertynumber)%npoints,2))
        this%spring12(i)%c%data_points(:,:) = spring12properties(this%spring12(i)%c%propertynumber)%data_points
      end do

    end if
    close(unit = 3100)
    
 ! coupling element
    this%nnodes_matrix12 = 0
    this%number_matrix12 = 0
    print *, "Parsing matrices for 12 DOF ..." ! print outputs only in the terminal
    open(unit = 3301, file = 'matrix12input.txt', status = 'old', iostat = io_error) 
    !old: The file already exists (nonexistence is an error)
    !iostat: zero if no error occurs, otherwise some positive number
   
    if (io_error == 0) then  ! if the file exists, then read the information.

    ! Reading the first three commented lines of the file
      do i = 1, linejump ! linejump=3 ; reads the first three lines of the file 
        read(3301, *)    ! read(unit number, *format type not specified), *: size of the string not defined in advance
      end do     
      
    ! Reading the 4th line and filling in the integer variables
      read(3301,*) this%nnodes_matrix12, this%number_matrix12 
      
      do i = 1, linejump
        read(3301,*)
      end do 

    ! Reading all nodenumber/matrix property combinations (combination = object) 
      allocate(this%matrix12%matrix12_object(this%nnodes_matrix12))
      do i = 1, this%nnodes_matrix12
        read(3301,*) this%matrix12%matrix12_object(i)%matrix12_strain_flag, this%matrix12%matrix12_object(i)%node1, this%matrix12%matrix12_object(i)%node2, this%matrix12%matrix12_object(i)%nproperty
      end do

    ! Three commented lines
      do i = 1, linejump
          read(3301,*)
      end do   
    
      allocate(this%matrix12%matrix12_property(this%number_matrix12))
    ! reading the matrix properties
      do i = 1, this%number_matrix12
        read(3301,*) this%matrix12%matrix12_property(i)%matrix12_matrix
        read(3301,*) this%matrix12%matrix12_property(i)%matrix12_m
        read(3301,*) this%matrix12%matrix12_property(i)%alpha_s, this%matrix12%matrix12_property(i)%alpha_v
        read(3301,*) this%matrix12%matrix12_property(i)%alpha_obj
        read(3301,*)
        read(3301,*)
        read(3301,*)
      end do 
   
    end if
    close(unit = 3301)
    
!< reading damping_12
    this%ndamping12 = 0
    this%ndamping12properties = 0
    print *, "Parsing springs for 12 DOF ..."
    open(unit = 3200, file = 'damping12input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file exists, then read the information.

      do i = 1, linejump
        read(3200, *)
      end do

      read(3200, *) this%ndamping12, this%ndamping12properties

      do i = 1, linejump
        read(3200, *)
      end do

      allocate(this%damping12(this%ndamping12))
      do i = 1, this%ndamping12
        read(3200, *) damping12_sort, damping12node, damping12property, damping12dir(:)

        select case (trim(adjustl(damping12_sort)))
          case ('angularvelocity_local')
            allocate(this%damping12(i)%c, source = damping12_angularvelocity_local() )

          case ('angularvelocity_global')
            allocate(this%damping12(i)%c, source = damping12_angularvelocity_global() )

          case ('translatory_local')
            allocate(this%damping12(i)%c, source = damping12_translatory_local() )

          case default
              print *, "ERROR: didn't recognize damping12: ", damping12_sort
              stop
          end select

        allocate(this%damping12(i)%c%sort, source = trim(adjustl(damping12_sort)))
        this%damping12(i)%c%node                  = damping12node
        this%damping12(i)%c%propertynumber        = damping12property
        this%damping12(i)%c%dir(:)                = damping12dir(:)

      end do

      allocate(damping12properties(this%ndamping12properties))
      do i = 1,this%ndamping12properties
        do j = 1, linejump
          read(3200, *)
        end do

        read(3200, *) damping12properties(i)%npoints
        allocate(damping12properties(i)%data_points(damping12properties(i)%npoints,2))
        do j = 1,damping12properties(i)%npoints
          read(3200, *) damping12properties(i)%data_points(j,:)
        end do
      end do

      do i = 1,this%ndamping12
        allocate(this%damping12(i)%c%data_points(damping12properties(this%damping12(i)%c%propertynumber)%npoints,2))
        this%damping12(i)%c%data_points(:,:) = damping12properties(this%damping12(i)%c%propertynumber)%data_points
      end do

    end if
    close(unit = 3200)

!< reading constraints6
    this%nconstraints6 = 0
    print *, "Parsing constraints for 6 DOF nodes..."
    open(unit = 4000, file = 'constraint6input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file constraint6input exists, then read the information.

       do i = 1, linejump
          read(4000, *)
       end do

       read(4000, *) this%nconstraints6
       allocate(this%constraints6(this%nconstraints6))

       do i = 1, linejump
          read(4000, *)
       end do

       do i = 1, this%nconstraints6
          read(4000, *) constraintsort, temp_constraint%nodes(:), temp_constraint%phi(:, 1), temp_constraint%phi(:, 2), temp_constraint%dir
          select case (constraintsort)
            case ('internal')
              allocate(this%constraints6(i)%c, source = constraint_6_internal() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint6input.txt! Definition of Node #2 is wrong!'
                temp_constraint%nodes(2) = 0
              end if

            case ('rotation_global')
              allocate(this%constraints6(i)%c, source = constraint_6_rotation_global() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint6input.txt! Definition of Node #2 is wrong!'
                temp_constraint%nodes(2) = 0
              end if

            case ('sphericalsupport')
              allocate(this%constraints6(i)%c, source = constraint_6_sphericalsupport() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint6input.txt! Definition of Node #2 is wrong!'
                temp_constraint%nodes(2) = 0
              end if

            case ('revolutesupport')
              allocate(this%constraints6(i)%c, source = constraint_6_revolutesupport() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint6input.txt! Definition of Node #2 is wrong!'
                temp_constraint%nodes(2) = 0
              end if

            case('sphericaljoint')
              allocate(this%constraints6(i)%c, source = constraint_6_sphericaljoint() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint6input.txt! Either Node #1 or Node #2 is zero!'
                stop
             end if

            case('revolutejoint')
              allocate(this%constraints6(i)%c, source = constraint_6_revolutejoint() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint6input.txt! Either Node #1 or Node #2 is zero!'
                stop
             end if

            case('inextensiblerevolutejoint')
              allocate(this%constraints6(i)%c, source = constraint_6_inextensiblerevolutejoint() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint6input.txt! Either Node #1 or Node #2 is zero!'
                stop
             end if

            case('layerconnection')
              allocate(this%constraints6(i)%c, source = constraint_6_layerconnection() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint6input.txt! Either Node #1 or Node #2 is zero!'
                stop
             end if

            case ('simplesupport')
              allocate(this%constraints6(i)%c, source = constraint_6_simplesupport() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint6input.txt! Definition of Node #2 is wrong!'
                temp_constraint%nodes(2) = 0
              end if

            case ('symmetry')
              allocate(this%constraints6(i)%c, source = constraint_6_symmetry() )

            case default
              print *, "ERROR: didn't recognize constraint: ", constraintsort
          end select

          this%constraints6(i)%c%nodes = temp_constraint%nodes
          this%constraints6(i)%c%phi   = temp_constraint%phi
          this%constraints6(i)%c%dir   = temp_constraint%dir
       end do
    end if
    close(unit = 4000)

!< reading constraints12
    this%nconstraints12 = 0
    print *, "Parsing constraints for 12 DOF nodes..."
    open(unit = 5000, file = 'constraint12input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file constraint12input exists, then read the information.
       do i = 1, linejump
          read(5000, *)
       end do

       read(5000, *) this%nconstraints12
       allocate(this%constraints12(this%nconstraints12))

       do i = 1, linejump
          read(5000, *)
       end do

       do i = 1, this%nconstraints12
          read(5000, *) constraintsort, temp_constraint%nodes(:), temp_constraint%phi(:, 1), temp_constraint%phi(:, 2), temp_constraint%dir
          select case (constraintsort)

            case ('internal')
              allocate(this%constraints12(i)%c, source = constraint_12_internal() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case ('internal_red')
              allocate(this%constraints12(i)%c, source = constraint_12_internal_red() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if
              
            case ('sphericalsupport')
              allocate(this%constraints12(i)%c, source = constraint_12_sphericalsupport() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case ('simplesupport_global')
              allocate(this%constraints12(i)%c, source = constraint_12_simplesupport_global() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case ('simplesupport_local')
              allocate(this%constraints12(i)%c, source = constraint_12_simplesupport_local() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case ('rigidsupport')
              allocate(this%constraints12(i)%c, source = constraint_12_rigidsupport() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case('sphericaljoint')
              allocate(this%constraints12(i)%c, source = constraint_12_sphericaljoint() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if

            case('revolutejoint')
              allocate(this%constraints12(i)%c, source = constraint_12_revolutejoint() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if

            case('rotorrevolutejoint')
              allocate(this%constraints12(i)%c, source = constraint_12_rotorrevolutejoint() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if

            case('revolutejoint_2')
              allocate(this%constraints12(i)%c, source = constraint_12_revolutejoint_2() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if
             
            case('revolutejoint_3')
              allocate(this%constraints12(i)%c, source = constraint_12_revolutejoint_3() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if

            case('rigidconnection')
              allocate(this%constraints12(i)%c, source = constraint_12_rigidconnection() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if

            case('masspointconnection')
              allocate(this%constraints12(i)%c, source = constraint_12_masspointconnection() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if

            case ('angularvelocity_global')
              allocate(this%constraints12(i)%c, source = constraint_12_angularvelocity_global() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case ('angularvelocity_local')
              allocate(this%constraints12(i)%c, source = constraint_12_angularvelocity_local() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case ('rotation_local')
              allocate(this%constraints12(i)%c, source = constraint_12_rotation_local() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if
              
            case ('relative_rotation_local')
              allocate(this%constraints12(i)%c, source = constraint_12_relative_rotation_local() )
              if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint12input.txt! Either Node #1 or Node #2 is zero in !', i
            end if

            case ('rotation_global')
              allocate(this%constraints12(i)%c, source = constraint_12_rotation_global() )
              if (temp_constraint%nodes(2) .ne. 0) then
                print*,'Warning: Check constraint12input.txt! Definition of Node #2 is wrong in !', i
                temp_constraint%nodes(2) = 0
              end if

            case default
              print *, "ERROR: didn't recognize constraint: ", constraintsort
          end select

          this%constraints12(i)%c%nodes           = temp_constraint%nodes
          this%constraints12(i)%c%phi             = temp_constraint%phi
          this%constraints12(i)%c%dir             = temp_constraint%dir

       end do
    end if
    close(unit = 5000)

!< reading constraint6to12
    this%nconstraints6to12 = 0
    print *, "Parsing constraints for 6 DOF and 12 DOF nodes..."
    open(unit = 5500, file = 'constraint6to12input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file constraint12input exists, then read the information.
       do i = 1, linejump
          read(5500, *)
       end do
       read(5500, *) this%nconstraints6to12
       allocate(this%constraints6to12(this%nconstraints6to12))

       do i = 1, linejump
          read(5500, *)
       end do

       do i = 1, this%nconstraints6to12
          read(5500, *) constraintsort, temp_constraint%nodes(:)
          select case (constraintsort)
            case ('softtransition')
              allocate(this%constraints6to12(i)%c, source = constraint_6to12_softtransition() )
             if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint6to12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
             end if

            case ('rigidtransition')
              allocate(this%constraints6to12(i)%c, source = constraint_6to12_rigidtransition() )
              if (temp_constraint%nodes(1) .eq. 0 .or. temp_constraint%nodes(2) .eq. 0) then
                print*,'Warning: Check joint connection in constraint6to12input.txt! Either Node #1 or Node #2 is zero in !', i
                stop
              end if

            case default
              print *, "ERROR: didn't recognize constraint: ", constraintsort
          end select
          this%constraints6to12(i)%c%nodes = temp_constraint%nodes
       end do
    end if
    close(unit = 5500)
    
    this%nboundaries12 = 0
    print *, "Parsing boundary condition amplitudes for 12 DOF constraints..."
    open(unit = 10000, file = 'boundary12input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file boundary12input exists, then read the information.
       do i = 1, linejump
          read(10000, *)
       end do

       read(10000, *) this%nboundaries12

       allocate(this%boundaries12(this%nboundaries12))
       allocate(this%boundaryamplitudes12(this%nboundaries12))

       do i = 1, linejump
          read(10000, *)
       end do

       do i = 1, this%nboundaries12
          read(10000, *) loadsort, this%boundaryamplitudes12(i)%intensity, this%boundaryamplitudes12(i)%duration, this%boundaries12(i)%constraint12, filename_boundary12
          allocate(this%boundaryamplitudes12(i)%sort, source = trim(loadsort))

          if (trim(filename_boundary12) .ne. 'nofile') then
            allocate(this%boundaryamplitudes12(i)%filename, source = trim(filename_boundary12))
          end if

       end do
    end if
    close(unit = 10000)

    this%nboundaries6 = 0
    open(unit = 10000, file = 'boundary6input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file boundary12input exists, then read the information.
       do i = 1, linejump
          read(10000, *)
       end do

       read(10000, *) this%nboundaries6
       allocate(this%boundaries6(this%nboundaries6))
       allocate(this%boundaryamplitudes6(this%nboundaries6))

       do i = 1, linejump
          read(10000, *)
       end do

       do i = 1, this%nboundaries6
          read(10000, *) loadsort, this%boundaryamplitudes6(i)%intensity, this%boundaryamplitudes6(i)%duration, this%boundaries6(i)%constraint6, filename_boundary6
          allocate(this%boundaryamplitudes6(i)%sort, source = trim(loadsort))

          if (trim(adjustl(filename_boundary6)) .ne. 'nofile') then
            allocate(this%boundaryamplitudes6(i)%filename, source = trim(adjustl(filename_boundary6)))
          end if

       end do
    end if
    close(unit = 10000)

!< reading loads6
    this%nloads6 = 0
    print *, "Parsing loads for 6 DOF constraints..."
    open(unit = 6000, file = 'load6input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file load6input exists, then read the information.
       do i = 1, linejump
          read(6000, *)
       end do
       read(6000, *) this%nloads6
       allocate(this%loads6(this%nloads6))
       allocate(this%loadamplitudes6(this%nloads6))
       do i = 1, linejump
          read(6000, *)
       end do
       do i = 1, this%nloads6
          read(6000, *) loadsort, this%loadamplitudes6(i)%intensity, this%loadamplitudes6(i)%duration, this%loads6(i)%node, this%loads6(i)%spatial(:), this%loads6(i)%material
          allocate(this%loadamplitudes6(i)%sort, source = trim(loadsort))
       end do
    end if
    close(unit = 6000)

!< reading loads12
    this%nloads12 = 0
    print *, "Parsing loads for 12 DOF constraints..."
    open(unit = 7000, file = 'load12input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file load6input exists, then read the information.
       do i = 1, linejump
          read(7000, *)
       end do
       read(7000, *) this%nloads12
       allocate(this%loads12(this%nloads12))
       allocate(this%loadamplitudes12(this%nloads12))
       do i = 1, linejump
          read(7000, *)
       end do
       do i = 1, this%nloads12
          read(7000, *) loadsort, this%loadamplitudes12(i)%intensity, this%loadamplitudes12(i)%duration, this%loads12(i)%node, this%loads12(i)%spatial(:), this%loads12(i)%material(:)
          allocate(this%loadamplitudes12(i)%sort, source = trim(loadsort))
       end do
    end if
    close(unit = 7000)

!< reading pmass6
    this%npmass6 = 0
    print *, "Parsing point masses for 6 DOF constraints..."
    open(unit = 8000, file = 'pointmass6input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file pmass6input exists, then read the information.
       do i = 1, linejump
          read(8000, *)
       end do
       read(8000, *) this%npmass6
       allocate(this%pmass6(this%npmass6))
       do i = 1, linejump
          read(8000, *)
       end do
       do i = 1, this%npmass6
          read(8000, *)  this%pmass6(i)%node, this%pmass6(i)%mass
       end do
    end if
    close(unit = 8000)

!< reading pmass12
    this%npmass12 = 0
    print *, "Parsing point masses for 12 DOF constraints..."
    open(unit = 9000, file = 'pointmass12input.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then  ! if the file pmass12input exists, then read the information.
       do i = 1, linejump
          read(9000, *)
       end do
       read(9000, *) this%npmass12
       allocate(this%pmass12(this%npmass12))
       do i = 1, linejump
          read(9000, *)
       end do
       do i = 1, this%npmass12
          read(9000, *)  this%pmass12(i)%node, this%pmass12(i)%mass
       end do
    end if
    close(unit = 9000)

!< reading simulation parameters
    print *, "Parsing simulation input file..."
    open(unit = 100, file = 'simulationinput_structure.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then ! if the file simulationinput exists, then read the information.
      ! reading filename
      do i = 1, linejump
          read(100, *)
      end do
      read(100, *) strfilename, stroldfilename
      allocate(this%output_filename, source =  trim(adjustl(strfilename)))
      if (trim(adjustl(stroldfilename)) .ne. 'none' .and. trim(adjustl(stroldfilename)) .ne. '!!') then
        this%boolean_oldsim = .TRUE.
        allocate(this%oldsimfilename, source =  trim(adjustl(stroldfilename)))
      end if
      if (trim(adjustl(stroldfilename)) .eq. '!!') backspace(100)
      
      ! reading simulation steps      
      do i = 1, linejump
          read(100, *)
      end do
      read(100, *) nStepCount
      do i = 1, linejump
          read(100, *)
      end do
      allocate(this%settings(nStepCount))
      do i = 1, nStepCount
         read(100, *) stepSort
         if (stepsort == 'dynamic' .or. stepsort == 'static') then
           read(100, *) this%settings(i)%totalt, this%settings(i)%deltat, this%settings(i)%tolerance, this%settings(i)%iterationlimit, this%settings(i)%gravityflag
           read(100, *) this%settings(i)%output_flag 
           allocate(this%settings(i)%simutype , source = trim(stepsort))
           do j = 1, linejump
              read(100, *)
           end do
         elseif (stepsort == 'static_arc') then
           read(100, *) this%settings(i)%totalt, this%settings(i)%deltat, this%settings(i)%tolerance, this%settings(i)%iterationlimit, this%settings(i)%arc_len_method_flag, this%settings(i)%number_desired_iteration
           read(100, *) this%settings(i)%output_flag
           allocate(this%settings(i)%simutype , source = trim(stepsort))
           do j = 1, linejump
              read(100, *)
           end do
         elseif (stepsort == 'modal' .or. stepsort == 'buckling') then
           read(100, *) this%settings(i)%nEF, this%settings(i)%tolerance, this%settings(i)%emin, this%settings(i)%emax, this%settings(i)%sparse_flag
           read(100, *) this%settings(i)%output_flag
           allocate(this%settings(i)%simutype , source = trim(stepsort))
           do j = 1, linejump
              read(100, *)
           end do
         elseif (stepsort == 'invariants') then
           read(100, *) this%settings(i)%gravityflag
           read(100, *) this%settings(i)%output_flag
           allocate(this%settings(i)%simutype , source = trim(stepsort))
           do j = 1, linejump
              read(100, *)
           end do
         else
           print*,'Warning: Simulation setting unknown, check simulationinput_structure.txt...'
           call exit(-1)
         end if
         this%settings(i)%step = i
      end do
      read(100, *) this%gravity(:)
      read(100, *, iostat=io)  this%settings(1)%optional_adaptive_deltat
        if (io .ne. 0) then 
            this%settings(1)%optional_adaptive_deltat = 0
        end if
    else
      print*,'Warning: no simulation file for DeSiO-Structure available...'
      this%boolean_Abort = .TRUE.
    end if
    close(100)

!< reading postprocessing parameter
    print *, "Parsing postprocessing input file..."
    open(unit = 110, file = 'postprocessinput_structure.txt', status = 'old', iostat = io_error)
    if (io_error == 0) then ! if the file postprocessinput_structure.txt exists, then read the information.
      this%postproc%boolean_q  = .TRUE.
      this%postproc%boolean_v  = .TRUE.
      this%postproc%boolean_lambda = .TRUE.
      this%postproc%boolean_e  = .FALSE.
      this%postproc%boolean_m = .FALSE.
      do i = 1, linejump
          read(110, *)
      end do
      read(110, *) npostproc
      do i = 1, linejump
          read(110, *)
      end do
      allocate(this%postproc(npostproc))
      do i = 1, npostproc
         read(110, *, iostat = io_error )  strtype
         if (io_error .ne. 0) then
            print *, "Error in postprocessing file!"
            call exit(-1)
         end if
         allocate(this%postproc(i)%strtype, source = trim(adjustl(strtype)))
         !allocate(this%postproc(i)%id, source = postprocid(1:count(postprocid .ne. 0))) 
      end do

    else
      print *, "Warning: no postprocessing input file found!"
    end if
    close(110)
    
    this%boolean_model_read = .TRUE.
    return
end subroutine readinginputs

!< Subroutine for reading old simulation files and continuing
  subroutine read_oldsim(this)
    implicit none
    class(model_structure), intent(inout) :: this
    integer :: io_error, i, j, nsteps
    real(kind = 8) :: time, temp_t, temp_lmomentum(3), temp_amomentum(3), temp_kenergy, temp_penergy, temp_tenergy, temp_val(3)
    real(kind = 8), allocatable :: temp_q(:), temp_v(:), temp_lambda(:)
    character(len = 1024) :: char_t
    logical :: bool_error = .FALSE.
    
    allocate(temp_q(size(this%q_t)))
    allocate(temp_v(size(this%qdot_t)))
    allocate(temp_lambda(size(this%lambda_t)))
    
    temp_q      = 0.0d0
    temp_v      = 0.0d0
    temp_lambda = 0.0d0
    
    !< reading result files from previous simulation
      open(unit = 5000, file = this%oldsimfilename // '_t.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error) go to 10
      
      !< determine number of rows in result files
      nsteps = 0
      read(5000, *, iostat = io_error)
      do 
        if (io_error .ne. 0) exit
        read(5000, *, iostat = io_error)
        nsteps = nsteps + 1
      end do
      close(unit = 5000)
      
      open(unit = 5000, file = this%oldsimfilename // '_t.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0) then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_t.dres'
        go to 10
      end if
        
      open(unit = 5001, file = this%oldsimfilename // '_q.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0) then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_q.dres'
        go to 10
      end if
      
      open(unit = 5002, file = this%oldsimfilename // '_v.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0) then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_v.dres'
        go to 10
      end if
      
      open(unit = 5003, file = this%oldsimfilename // '_lambda.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0) then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_lambda.dres'
        go to 10
      end if
      
      !open(unit = 5004, file = this%oldsimfilename // '_e.dres', status = 'old', iostat = io_error)
      !if (io_error .ne. 0) bool_error = .TRUE.
      !if (bool_error) go to 10
      
      print*, '... all structure simulation files read correctly'
      
      ! jumping to last row-1 in previos result files 
      do i = 1, nsteps-1
        read(5000, *, iostat = io_error)
        if (io_error .ne. 0) bool_error = .TRUE.
        if (bool_error) go to 10
        
        read(5001, *, iostat = io_error)
        if (io_error .ne. 0) bool_error = .TRUE.
        if (bool_error) go to 10
        
        read(5002, *, iostat = io_error)
        if (io_error .ne. 0) bool_error = .TRUE.
        if (bool_error) go to 10
        
        read(5003, *, iostat = io_error)
        if (io_error .ne. 0) bool_error = .TRUE.
        if (bool_error) go to 10
        
        !read(5004, *, iostat = io_error)
        !if (io_error .ne. 0) bool_error = .TRUE.
        !if (bool_error) go to 10
      end do
      
      ! reading last row of previos result files
      read(5000, *, iostat = io_error) temp_t, time
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error) go to 10
      
      read(5001, *, iostat = io_error) temp_q
      if ( io_error .ne. 0 ) bool_error = .TRUE.
      if (bool_error) go to 10
      
      read(5002, *, iostat = io_error) temp_v
      if ( io_error .ne. 0 ) bool_error = .TRUE.
      if (bool_error) go to 10
            
      read(5003, *, iostat = io_error) temp_lambda
      if ( io_error .ne. 0 ) bool_error = .TRUE.
      if (bool_error) go to 10

      !read(5004, *, iostat = io_error) temp_lmomentum(1:3), temp_amomentum(1:3), temp_kenergy, temp_penergy, temp_tenergy, temp_val(1:3)
      !if ( io_error .ne. 0 ) bool_error = .TRUE.
      !if (bool_error) go to 10

      !< assign surface results
      this%q_t(:)      = temp_q(1:size(this%q_t))
      this%qdot_t(:)   = temp_v(1:size(this%qdot_t))
      this%lambda_t(:) = temp_lambda(1:size(this%lambda_t))
      this%time_t      = temp_t
      this%kenergy_t   = 0.0d0
      this%penergy_t   = 0.0d0
      this%tenergy_t   = 0.0d0
      this%lmomentum_t = 0.0d0
      this%amomentum_t = 0.0d0
      
      !< if error, then continue from here
10    continue
      
      close(unit = 5000)
      close(unit = 5001)
      close(unit = 5002)
      close(unit = 5003)
      !close(unit = 5004)
      
      if (bool_error) then
        print*, 'Warning in structural simulation! No result or error in files found! Simulation restarts!'
        this%boolean_oldsim = .FALSE.
      end if
      
    return
  end subroutine read_oldsim

!< Subroutine for sorting a vector in ascending order using the linear sorting algorithm
  subroutine linSorting_int(inzToSort,arrayToSort)
    implicit none

    logical :: bConverged

    integer :: i, j, k, r1
    integer, allocatable, dimension(:), intent(inout) :: inzToSort
    integer :: inzToCheck, inzToCheck_t

    integer, allocatable, dimension(:), intent(inout) :: arrayToSort
    integer :: valToCheck, valToCheck_t

    allocate(inzToSort(size(arrayToSort)))
    inzToSort(1:size(arrayToSort)) = (/ (r1, r1 = 1,size(arrayToSort)) /)

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

!< Subroutine to calculate the null space of a matrix using the Singular Value Decomposition
  subroutine get_nullspace_mod(nullspace_matrix,matrix,nq,nc)
    real(kind = 8), allocatable :: SV(:,:), VT(:,:), U(:,:), work(:), matrix_temp(:,:)
    real(kind = 8), allocatable, intent(inout) :: nullspace_matrix(:,:)
    real(kind = 8), allocatable :: tr_nullspace_matrix(:,:)
    real(kind = 8), intent(inout) :: matrix(:,:)
    integer :: info, nq, nc, r1, lwork
    integer, allocatable :: indices_vec(:)

    lwork = (max(3*min(nc,nq) + max(nc,nq),5*min(nc,nq)))
    
    allocate(work(lwork))
    allocate(SV(nc,nq))
    allocate(VT(nq,nq))
    allocate(U(nc,nc))
    allocate(matrix_temp(nc,nq))
    
    matrix_temp = matrix
    SV = 0.0d0
    VT = 0.0d0
    U  = 0.0d0

    call DGESVD( 'N', 'A',nc, nq, matrix_temp, nc, SV, U, 1, VT, nq, work, lwork, info)

! getting nullspace of matrix
    if (allocated(nullspace_matrix)) then
      deallocate(nullspace_matrix)
    end if

    allocate(tr_nullspace_matrix(nq-nc,nq))
    allocate(nullspace_matrix(nq,nq-nc))
    allocate(indices_vec(nq-nc))

    indices_vec(1:nq-nc)     = (/ (r1, r1 = nc+1,nq) /)
    tr_nullspace_matrix      = 0.0d0
    tr_nullspace_matrix(:,:) = VT(indices_vec,:)
    nullspace_matrix(:,:)    = transpose(tr_nullspace_matrix)

    deallocate(tr_nullspace_matrix)
    deallocate(work)
    deallocate(SV)
    deallocate(VT)
    deallocate(U)
    deallocate(matrix_temp)
  end subroutine get_nullspace_mod

!< Subroutine to append a value to a matrix
  subroutine reallocate_IArray(array,n_i,m_j,val)
  implicit none

  integer, intent(in) :: n_i, m_j
  integer :: n_a, m_a, n, m
  integer :: i, j
  integer, intent(in) :: val
  integer, allocatable, intent(inout) :: array(:,:)
  integer, allocatable :: temp_array(:,:)

  if (.not. allocated(array)) then
    allocate(array(n_i,m_j))
    array(n_i,m_j) = val
    return
  end if

  n_a = size(array,1)
  m_a = size(array,2)

  allocate(temp_array(n_a,m_a))
  temp_array(:,:) = array(:,:)

  deallocate(array)

  n = n_a
  if (n_i>n_a) then
    n = n_i
  end if

  m = m_a
  if (m_j>m_a) then
    m = m_j
  end if

  allocate(array(n,m))
  array(:,:) = 0

  do i = 1,n_a
      do j = 1,m_a
        array(i,j) = temp_array(i,j)
      end do
  end do
  array(n_i,m_j) = val

end subroutine reallocate_IArray

  subroutine extract_csr_Format_From_Sparsematrix(sparsematrix, csr_values, csr_columns,csr_rowIndices, csr_compr_source, indicesRow,indicesCol)
!< output indices (rows, columns) and values are reordered new
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

        vec_range(1:size_vec_range) = (/ (r1, r1 = csr_RowIndices_inp(indicesRow(i)),csr_RowIndices_inp(indicesRow(i)+1)-1) /)
        values_n  = csr_Values_inp(vec_range)
        columns_n = csr_Columns_inp(vec_range)

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

!< Subroutine for writing character line from file to an array/vector
  subroutine read_line_to_array(file_unit, vec, var_type, ireason)
    implicit none
    
    integer, intent(in):: file_unit
    integer, intent(inout):: ireason
    integer :: k, var_type
    
    character(len = 1) :: strline
    character(len=:), allocatable :: temp_strLine
    
    integer :: val
    integer, allocatable, intent(inout) :: vec(:)
    
    logical :: boolean_valfound = .FALSE.

    temp_strLine = ''
    if (allocated(vec)) deallocate(vec)
    k = 0
    do
      read(file_unit,'(a)',advance='NO',iostat=ireason) strline
      
      if (k .eq. 0) then
        if (ireason .ne. 0) then
            ireason = -1
            exit
        end if
      end if
      
      if (strline .ne. ' ') then
        boolean_valfound = .TRUE.
        temp_strLine = trim(adjustl(temp_strLine)) // trim(adjustl(strline))
      else
        if (boolean_valfound) then
          if (var_type == 0) then
            read(temp_strLine, '(f30.15)' )  val
            call append_int_vec(vec, val) 
          end if
          
        end if
        temp_strLine = ''
        boolean_valfound = .FALSE.
      end if
      
      if (IS_IOSTAT_END(ireason)) then
        ireason = 0
        exit
      end if

      k = k + 1
    end do
          
    return
  end subroutine read_line_to_array
  
end module class_model_structure
