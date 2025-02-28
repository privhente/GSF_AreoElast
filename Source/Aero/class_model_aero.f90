module class_model_aero

  use my_constants_aero, only: pi
  use my_types_aero, only: tangent_matrices3_3x3
  use class_bounded_vortex_sheet
  use class_unbounded_vortex_sheet
  use class_inflow_aero
  use my_vlm_functions, only: append_real_vec, read_line_to_array, find_segments, append_int_vec
  use my_FileIO_functions
  use sparse_function, only: single_sparse_matrix, create_rowcol_format, sparse_initialize
  implicit none

  type :: model_aero

    type(bounded_vortex_sheet), allocatable :: surfaces(:)
    type(unbounded_vortex_sheet), allocatable :: wakes(:)
    type(air_flow) :: airflow
    type(single_sparse_matrix) :: sparse_Tcp1, sparse_Tcp3, sparse_dfae

    !< arrays for model aerodynamic simulation
    integer :: nsurfaces
    integer :: nwakes
    integer :: nwakesproperties
    real(kind = 8), allocatable :: amatrix(:, :), inv_amatrix(:,:)
    real(kind = 8), allocatable :: circulations(:), oldcirculations(:)
    real(kind = 8), allocatable :: rhs(:)
    real(kind = 8), allocatable :: dAG(:,:), drhs(:,:), dG(:,:)
    real(kind = 8), allocatable :: dfae(:,:), dfae_ring(:,:)
    real(kind = 8), allocatable :: qs_t(:), vs_t(:)
    real(kind = 8), allocatable :: deltaps(:), fae(:), faen(:), areanormal(:), fae_x(:), fae_w(:), fae_xn(:), fae_wn(:)
    real(kind = 8), allocatable :: deltav(:), vmean(:), vs(:), vsv(:), vwv(:)

    integer :: tnrings
    integer :: tnnodes, tnnodew
    integer :: tncoordinates, tnvelocities

    integer :: nsteps = 0
    integer :: tcounter = 0
    real(kind = 8) :: totalt, deltat, time
    real(kind = 8) :: cutoff
    real(kind = 8) :: delta_sim_time = 0.0d0
    real(kind = 8), allocatable, dimension(:,:) :: arrqs_t, arrvs_t
    character(:), allocatable :: output_filename, oldsimfilename, strSimuType

    !< arrays for output
    real(kind = 8), allocatable :: circulations_nodal(:), deltaps_nodal(:), areanormal_nodal(:), vmean_nodal(:)

    !< array for sensitivity analysis
    real(kind = 8), allocatable :: arrgradqs_t(:,:), arrgradvs_t(:,:), gradqsvs_t(:)
    real(kind = 8), allocatable :: deltax_t(:)

    !< more variable for solver
    integer, allocatable :: ipiv(:)
    real(kind = 8), allocatable :: work(:)

    logical :: bool_const_a_matrix = .TRUE.
    logical :: boolean_abort = .FALSE.
    logical :: boolean_update_simulation = .TRUE.
    logical :: boolean_model_read = .FALSE.
    logical :: boolean_oldsim = .FALSE.
    logical :: boolean_actualize_amatrix = .TRUE.
    logical :: boolean_unsteady_term = .TRUE.
    logical :: boolean_abort_airflow = .FALSE.
    logical :: boolean_lin_wake_row = .True.
  contains

     procedure :: constructor
     procedure :: actualize
     procedure :: evaluateamatrix
     procedure :: evaluaterhs
     procedure :: update_wake
     procedure :: evaluateamplitude
     procedure :: read_kinematic
     procedure :: output_open
     procedure :: output_write
     procedure :: output_close
     procedure :: wake_velocity
     procedure :: wake_init
     procedure :: solver
     procedure :: modelfae
     procedure :: modeldfae
     procedure :: read_oldsim
     procedure :: readinginputs
     procedure :: read_oldwake
     procedure :: solver_num_sensitivity
     procedure :: wake_update_first_row
     procedure :: linearized_non_penetration_condition

  end type model_aero

  contains

!< Subroutine for construct the model
  subroutine constructor(this)
    implicit none

    class(model_aero), intent(inout) :: this
    integer :: i, j, k, l, j3w, j3cs, j3qs, j3ps, j3a, j3b, j3c, j3d, j1ns
    integer :: j3cw, j3qw, j3pw, j1nw
    integer, dimension(2) :: segmentconnectivitydir, segmentconnectivityinv
    integer :: indicesp_ae_global(3), indicesc_ae_global(1)
    integer, dimension(1) :: arr_j1, arr_j2, arr_j3, arr_j4
    real(kind = 8) :: arr_val(1,1)
    integer, allocatable :: arrayaux1(:, :), arrayaux2(:, :), arrayaux3(:, :), temp_rings(:)
    integer :: r1
    integer :: convert_job_coo_csr(8), info, nNonZeros
    integer :: enode, aa, n_wake_edge_segments, ID_first_row_wake_edge, node1, node2
    integer, allocatable :: temp_vec(:)
    
    integer :: segment1_number, segment2_number, segment3_number, segment4_number
    
    convert_job_coo_csr    = 0 ! the matrix in the CSR format is converted to the coordinate format;
    convert_job_coo_csr(1) = 2 ! the matrix in the coordinate format is converted to the CSR format, and the column indices in CSR representation are sorted in the increasing order within each row.
    convert_job_coo_csr(2) = 1 ! one-based indexing for the matrix in CSR format is used.
    convert_job_coo_csr(3) = 1 ! one-based indexing for the matrix in coordinate format is used.
    convert_job_coo_csr(6) = 0 ! all arrays acsr, ja, ia are filled in for the output storage.

    this%nsurfaces = 0
    if (allocated(this%surfaces)) then
        this%nsurfaces = size(this%surfaces)
    end if

    this%nwakes = 0
    if (allocated(this%wakes)) then
        this%nwakes = size(this%wakes)
    end if

    if (this%nsurfaces .eq. 0) then
      this%boolean_abort = .TRUE.
      return
    end if

    this%tnrings       = 0
    this%tnnodes       = 0
    this%tncoordinates = 0
    this%tnvelocities  = 0
    this%tnnodew       = 0

!< initalizing and allocating surfaces
    j3cs = 0
    j3qs = 0
    j3ps = 0
    j1ns = 0
    print*, 'Allocating bounded vortex elements'
    do i = 1, this%nsurfaces
      !< allocating segments
      this%surfaces(i)%nsegments = this%surfaces(i)%nnodes1*(this%surfaces(i)%nnodes2-1)+(this%surfaces(i)%nnodes1-1)*this%surfaces(i)%nnodes2
      allocate(this%surfaces(i)%segments(this%surfaces(i)%nsegments))

      !< allocating segment cirucaltion
      allocate(this%surfaces(i)%segmentscirculations(this%surfaces(i)%nsegments))
      allocate(this%surfaces(i)%oldsegmentscirculations(this%surfaces(i)%nsegments))

      allocate(arrayaux1(this%surfaces(i)%nsegments, 2))
      allocate(arrayaux2(this%surfaces(i)%nsegments, 2))
      allocate(arrayaux3(this%surfaces(i)%nrings   , 4))

      !< define connectivity and adjencies for segments and define rings of the grid
      call find_segments(this%surfaces(i)%nnodes1, this%surfaces(i)%nnodes2, arrayaux1, arrayaux2, arrayaux3)

      !print*, '... Getting segments: seg. nr(1); node conn.(2-3); adj. ring(4-5)'
      do j = 1, this%surfaces(i)%nsegments
        this%surfaces(i)%segments(j)%connectivity  = arrayaux1(j, :)
        this%surfaces(i)%segments(j)%adjacency     = arrayaux2(j, :)
        this%surfaces(i)%segments(j)%indicesq(1:3) = (/ (r1, r1 = 3*(this%surfaces(i)%segments(j)%connectivity(1)-1)+1, 3*(this%surfaces(i)%segments(j)%connectivity(1)-1)+3) /)
        this%surfaces(i)%segments(j)%indicesq(4:6) = (/ (r1, r1 = 3*(this%surfaces(i)%segments(j)%connectivity(2)-1)+1, 3*(this%surfaces(i)%segments(j)%connectivity(2)-1)+3) /)
        !print*, j, this%surfaces(i)%segments(j)%connectivity, this%surfaces(i)%segments(j)%adjacency
      end do

      !print*, '... Getting circulation rings: ring nr.(1); adj. seg.(2-5)'
      do j = 1, this%surfaces(i)%nrings
        this%surfaces(i)%rings(j)%adjacency = arrayaux3(j, :)
        !print*, j, this%surfaces(i)%rings(j)%adjacency
      end do

      deallocate(arrayaux1)
      deallocate(arrayaux2)
      deallocate(arrayaux3)

      !< determine transormation matrix for transformation between ring circulation and (global) segment circulation
      call this%surfaces(i)%rings2segments_matrix()

      !< defining sections
      this%surfaces(i)%nsections = this%surfaces(i)%nnodes1-1
      allocate(this%surfaces(i)%sections(this%surfaces(i)%nsections))
      do j = 1, this%surfaces(i)%nsections
        this%surfaces(i)%sections(j)%nsegments = this%surfaces(i)%nnodes2
        allocate(this%surfaces(i)%sections(j)%segments(this%surfaces(i)%sections(j)%nsegments))
        do k = 1, this%surfaces(i)%sections(j)%nsegments
            this%surfaces(i)%sections(j)%segments(k) = this%surfaces(i)%nnodes2*(j-1)+k
        end do

        !< defining rings/circulation connectivity
        this%surfaces(i)%sections(j)%nrings = this%surfaces(i)%nnodes2-1
        allocate(this%surfaces(i)%sections(j)%rings(this%surfaces(i)%sections(j)%nrings))
        do k = 1, this%surfaces(i)%sections(j)%nrings
            this%surfaces(i)%sections(j)%rings(k) = j+(this%surfaces(i)%nnodes2-1)*(k-1)
        end do
      end do

      !< global surface variables
      this%tnrings       = this%tnrings       + this%surfaces(i)%nrings
      this%tnnodes       = this%tnnodes       + this%surfaces(i)%nnodes
      this%tncoordinates = this%tncoordinates + this%surfaces(i)%nnodes*3
      this%tnvelocities  = this%tnvelocities  + this%surfaces(i)%nnodes*3

      !< setting local indices for surface position, velocity, pressure and circulation
      do j = 1, this%surfaces(i)%nrings
        this%surfaces(i)%rings(j)%indicesp = (/ (r1, r1 = 3*(j-1)+1,3*(j-1)+3) /)
        this%surfaces(i)%rings(j)%indicesc = j

        j3a = 3*(this%surfaces(i)%rings(j)%connectivity(1)-1)
        j3b = 3*(this%surfaces(i)%rings(j)%connectivity(2)-1)
        j3c = 3*(this%surfaces(i)%rings(j)%connectivity(3)-1)
        j3d = 3*(this%surfaces(i)%rings(j)%connectivity(4)-1)

        !(/ (r1, r1 = 1,18) /)
        this%surfaces(i)%rings(j)%indicesq(1:3)   = (/ (r1, r1 = j3a+1,j3a+3) /)
        this%surfaces(i)%rings(j)%indicesq(4:6)   = (/ (r1, r1 = j3b+1,j3b+3) /)
        this%surfaces(i)%rings(j)%indicesq(7:9)   = (/ (r1, r1 = j3c+1,j3c+3) /)
        this%surfaces(i)%rings(j)%indicesq(10:12) = (/ (r1, r1 = j3d+1,j3d+3) /)

        !< indicesv and incicesq are related to the four nodes of ring: Mainly for determining control point ccordinates from nodal coordinates
        this%surfaces(i)%rings(j)%indicesv(1:12)  = this%surfaces(i)%rings(j)%indicesq(1:12)

        do k = 1,4
          this%surfaces(i)%rings(j)%nodes(k)%node_aero%localid  = this%surfaces(i)%rings(j)%connectivity(k)
          this%surfaces(i)%rings(j)%nodes(k)%node_aero%globalid = j1ns + this%surfaces(i)%rings(j)%connectivity(k)
        end do
      end do

      !< defining local indices for surface nodes
      do j = 1,this%surfaces(i)%nnodes
        this%surfaces(i)%nodes(j)%indicesq   = (/ (r1, r1 = 3*(j-1)+1,3*(j-1)+3) /)
        this%surfaces(i)%nodes(j)%indicesv   = (/ (r1, r1 = 3*(j-1)+1,3*(j-1)+3) /)
        this%surfaces(i)%nodes(j)%node_aero%localid  = j
        this%surfaces(i)%nodes(j)%node_aero%globalid = j1ns + j
      end do

      !< allocate ringcirculations, position and velocity of control points
      allocate(this%surfaces(i)%ringscirculations(this%surfaces(i)%nrings))
      allocate(this%surfaces(i)%oldringscirculations(this%surfaces(i)%nrings))
      allocate(this%surfaces(i)%q_t(3*this%surfaces(i)%nnodes))
      allocate(this%surfaces(i)%v_t(3*this%surfaces(i)%nnodes))

      this%surfaces(i)%ringscirculations    = 0.0d0
      this%surfaces(i)%oldringscirculations = 0.0d0
      this%surfaces(i)%q_t = 0.0d0
      this%surfaces(i)%v_t = 0.0d0

      !< global indices for circulation at rings, coordinates and  velocity at nodes
      allocate(this%surfaces(i)%indicesc(this%surfaces(i)%nrings))
      this%surfaces(i)%indicesc(1:this%surfaces(i)%nrings) = j3cs + (/ (r1, r1 = 1,this%surfaces(i)%nrings) /)

      allocate(this%surfaces(i)%indicesp(3*this%surfaces(i)%nrings))
      this%surfaces(i)%indicesp(1:3*this%surfaces(i)%nrings) = j3ps + (/ (r1, r1 = 1,3*this%surfaces(i)%nrings) /)

      allocate(this%surfaces(i)%indicesq(3*this%surfaces(i)%nnodes))
      this%surfaces(i)%indicesq(1:3*this%surfaces(i)%nnodes) = j3qs + (/ (r1, r1 = 1,3*this%surfaces(i)%nnodes) /)

      allocate(this%surfaces(i)%indicesv(3*this%surfaces(i)%nnodes))
      this%surfaces(i)%indicesv(1:3*this%surfaces(i)%nnodes) = j3qs + (/ (r1, r1 = 1,3*this%surfaces(i)%nnodes) /)

      j3cs = j3cs + this%surfaces(i)%nrings
      j3qs = j3qs + 3*this%surfaces(i)%nnodes
      j3ps = j3ps + 3*this%surfaces(i)%nrings
      j1ns = j1ns + this%surfaces(i)%nnodes
    end do

    ! initialize sparse matrices
    call sparse_initialize(this%sparse_Tcp1,   this%tnrings, this%tncoordinates)
    call sparse_initialize(this%sparse_Tcp3, 3*this%tnrings, this%tncoordinates)

    print*, 'Calculating sparse transformation matrix (nodes and cp)'
    !< set global sparse transformation matrix for tranfering cp to nodal values
    do i = 1, this%nsurfaces
      do j = 1, this%surfaces(i)%nrings
        allocate(this%surfaces(i)%rings(j)%avector(3*this%tnrings))

        indicesp_ae_global = this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp)
        call create_rowcol_format(this%sparse_Tcp3, indicesp_ae_global, this%surfaces(i)%indicesq(this%surfaces(i)%rings(j)%indicesq(1:3)),   0.25d0*i_33) !< node 1 of ring j
        call create_rowcol_format(this%sparse_Tcp3, indicesp_ae_global, this%surfaces(i)%indicesq(this%surfaces(i)%rings(j)%indicesq(4:6)),   0.25d0*i_33) !< node 2 of ring j
        call create_rowcol_format(this%sparse_Tcp3, indicesp_ae_global, this%surfaces(i)%indicesq(this%surfaces(i)%rings(j)%indicesq(7:9)),   0.25d0*i_33) !< node 3 of ring j
        call create_rowcol_format(this%sparse_Tcp3, indicesp_ae_global, this%surfaces(i)%indicesq(this%surfaces(i)%rings(j)%indicesq(10:12)), 0.25d0*i_33) !< node 4 of ring j

        this%surfaces(i)%rings(j)%Tcp(1:3,1:3)   = 0.25d0*i_33
        this%surfaces(i)%rings(j)%Tcp(1:3,4:6)   = 0.25d0*i_33
        this%surfaces(i)%rings(j)%Tcp(1:3,7:9)   = 0.25d0*i_33
        this%surfaces(i)%rings(j)%Tcp(1:3,10:12) = 0.25d0*i_33

        indicesc_ae_global = this%surfaces(i)%indicesc(this%surfaces(i)%rings(j)%indicesc)
        arr_j1  = this%surfaces(i)%rings(j)%nodes(1)%node_aero%globalID
        arr_j2  = this%surfaces(i)%rings(j)%nodes(2)%node_aero%globalID
        arr_j3  = this%surfaces(i)%rings(j)%nodes(3)%node_aero%globalID
        arr_j4  = this%surfaces(i)%rings(j)%nodes(4)%node_aero%globalID
        arr_val = 0.25d0
        call create_rowcol_format(this%sparse_Tcp1, indicesc_ae_global, arr_j1, arr_val)
        call create_rowcol_format(this%sparse_Tcp1, indicesc_ae_global, arr_j2, arr_val)
        call create_rowcol_format(this%sparse_Tcp1, indicesc_ae_global, arr_j3, arr_val)
        call create_rowcol_format(this%sparse_Tcp1, indicesc_ae_global, arr_j4, arr_val)
      end do
    end do

    !< convert transformation matrix Tcp3 from coo to csr sparse format
    allocate(this%sparse_Tcp3%csr_values( this%sparse_Tcp3%ncoords))
    allocate(this%sparse_Tcp3%csr_columns( this%sparse_Tcp3%ncoords))
    allocate(this%sparse_Tcp3%csr_rowIndices(3*this%tnrings+1))
    convert_job_coo_csr(5) =  this%sparse_Tcp3%ncoords ! maximum number of the non-zero elements allowed if job(1)=0.
    !< convert coo format to csr format: call mkl_scsrcoo(job, n, acsr, ja, ia, nnz, acoo, rowind, colind, info)
    call mkl_dcsrcoo(convert_job_coo_csr, 3*this%tnrings, this%sparse_Tcp3%csr_values, &
                  this%sparse_Tcp3%csr_columns, this%sparse_Tcp3%csr_rowIndices, this%sparse_Tcp3%ncoords, &
                  this%sparse_Tcp3%coo_values(1:this%sparse_Tcp3%ncoords), this%sparse_Tcp3%coo_rows(1:this%sparse_Tcp3%ncoords), this%sparse_Tcp3%coo_cols(1:this%sparse_Tcp3%ncoords), info)
    deallocate(this%sparse_Tcp3%coo_rows)
    deallocate(this%sparse_Tcp3%coo_cols)
    deallocate(this%sparse_Tcp3%coo_values)

    !< convert transformation matrix Tcp1 from coo to csr sparse format
    allocate(this%sparse_Tcp1%csr_values( this%sparse_Tcp1%ncoords))
    allocate(this%sparse_Tcp1%csr_columns( this%sparse_Tcp1%ncoords))
    allocate(this%sparse_Tcp1%csr_rowIndices(this%tnrings+1))
    convert_job_coo_csr(5) =  this%sparse_Tcp1%ncoords ! maximum number of the non-zero elements allowed if job(1)=0.
    !< convert coo format to csr format: call mkl_scsrcoo(job, n, acsr, ja, ia, nnz, acoo, rowind, colind, info)
    call mkl_dcsrcoo(convert_job_coo_csr, this%tnrings, this%sparse_Tcp1%csr_values, &
                  this%sparse_Tcp1%csr_columns, this%sparse_Tcp1%csr_rowIndices, this%sparse_Tcp1%ncoords, &
                  this%sparse_Tcp1%coo_values(1:this%sparse_Tcp1%ncoords), this%sparse_Tcp1%coo_rows(1:this%sparse_Tcp1%ncoords), this%sparse_Tcp1%coo_cols(1:this%sparse_Tcp1%ncoords), info)
    deallocate(this%sparse_Tcp1%coo_rows)
    deallocate(this%sparse_Tcp1%coo_cols)
    deallocate(this%sparse_Tcp1%coo_values)

!< identifying segments that correspond to an edge of the surface.
    do i = 1, this%nsurfaces
        do j = 1, this%surfaces(i)%nsegments
          if ((this%surfaces(i)%segments(j)%adjacency(1) == 0) .or. (this%surfaces(i)%segments(j)%adjacency(2) == 0)) then
              this%surfaces(i)%segments(j)%surfaceedge = 1
          end if
        end do
    end do

!< identifying segments of the bounded vortex sheet to the wake as wake edge
    do i = 1, this%nsurfaces
      do  j = 1, this%surfaces(i)%nsegments
        do k = 1, this%nwakes
            do l = 1, this%wakes(k)%nnodess-1
              segmentconnectivitydir = [this%wakes(k)%enodes(l), this%wakes(k)%enodes(l+1)]
              segmentconnectivityinv = [this%wakes(k)%enodes(l+1), this%wakes(k)%enodes(l)]
              if (all(this%surfaces(i)%segments(j)%connectivity == segmentconnectivitydir) .or. all(this%surfaces(i)%segments(j)%connectivity == segmentconnectivityinv)) then
                  this%surfaces(i)%segments(j)%wakeedge = 1
              end if
            end do
        end do
      end do
    end do

!< assigning edge types to ring segments
    do i = 1, this%nsurfaces
      do j = 1, this%surfaces(i)%nrings
        do k = 1, 4
          if (this%surfaces(i)%segments(this%surfaces(i)%rings(j)%adjacency(k))%surfaceedge == 1) then
            this%surfaces(i)%rings(j)%surfaceedge(k) = 1
          else
            this%surfaces(i)%rings(j)%surfaceedge(k) = 0
          end if
          if (this%surfaces(i)%segments(this%surfaces(i)%rings(j)%adjacency(k))%wakeedge == 1) then
            this%surfaces(i)%rings(j)%wakeedge(k) = 1
          else
            this%surfaces(i)%rings(j)%wakeedge(k) = 0
          end if
        end do
      end do
    end do

!< initializing and allocating wake surfaces
    j3w  = 0
    j3cw = 0
    j3qw = 0
    j3pw = 0
    j1nw = 0
    print*, 'Allocating unbounded vortex elements'
    do i = 1, this%nwakes
      !< assigning, if a surface is a lifting surface or not
      this%surfaces(this%wakes(i)%surface)%boolean_lifting_surface = .TRUE.

      ! determine number of steps to consider until wake is cut
      this%wakes(i)%nsteps       = this%nsteps
      this%wakes(i)%nstepstotal  = this%wakes(i)%nrowstotal
      this%wakes(i)%tcounter     = this%tcounter
      this%wakes(i)%nringst      = this%wakes(i)%nrowstotal
      this%wakes(i)%nnodest      = this%wakes(i)%nringst + 1
      this%wakes(i)%nnodes       = this%wakes(i)%nnodess*this%wakes(i)%nnodest
      this%wakes(i)%nrings       = this%wakes(i)%nringss*this%wakes(i)%nringst
      this%wakes(i)%ncoordinates = 3*this%wakes(i)%nnodes
      this%wakes(i)%nsegments    = this%wakes(i)%nringst*this%wakes(i)%nnodess + (this%wakes(i)%nnodess-1)*(this%wakes(i)%nringst+1)
      this%wakes(i)%nnodes1      = this%wakes(i)%nnodess
      this%wakes(i)%nnodes2      = this%wakes(i)%nnodest

      this%tnnodew = this%tnnodew + this%wakes(i)%nnodes
      allocate(this%wakes(i)%nodes(this%wakes(i)%nnodes))
      allocate(this%wakes(i)%rings(this%wakes(i)%nrings))
      allocate(this%wakes(i)%segments(this%wakes(i)%nsegments))

      allocate(this%wakes(i)%ringscirculations(this%wakes(i)%nrings))
      allocate(this%wakes(i)%temp_ringscirculations(this%wakes(i)%nrings))
      allocate(this%wakes(i)%segmentscirculations(this%wakes(i)%nsegments))
      allocate(this%wakes(i)%q_t(this%wakes(i)%ncoordinates))
      allocate(this%wakes(i)%q_t_first_row(3*this%wakes(i)%nnodess))
      allocate(this%wakes(i)%v_t(this%wakes(i)%ncoordinates))

      this%wakes(i)%ringscirculations(:) = 0.0d0
      this%wakes(i)%temp_ringscirculations(:) = 0.0d0
      this%wakes(i)%segmentscirculations(:) = 0.0d0
      this%wakes(i)%q_t(:) = 0.0d0
      this%wakes(i)%q_t_first_row = 0.0d0
      this%wakes(i)%v_t(:) = 0.0d0

      !< defining global indices for wake nodes
      do j = 1,this%wakes(i)%nnodes
        this%wakes(i)%nodes(j)%indicesq = (/ (r1, r1 = 3*(j-1)+1,3*(j-1)+3) /)
        this%wakes(i)%nodes(j)%indicesv = (/ (r1, r1 = 3*(j-1)+1,3*(j-1)+3) /)
      end do

      !< compute wake-node connectivity
      call this%wakes(i)%computeconnectivity()

      !< defining number of lines in the lattice
      allocate(arrayaux1(this%wakes(i)%nsegments, 2))
      allocate(arrayaux2(this%wakes(i)%nsegments, 2))
      allocate(arrayaux3(this%wakes(i)%nrings   , 4))

      !< define connectivity and adjencies for segments and define rings of the grid
      call find_segments(this%wakes(i)%nnodes1, this%wakes(i)%nnodes2, arrayaux1, arrayaux2, arrayaux3)

      !print*, '... Getting segments: seg. nr(1); node conn.(2-3); adj. ring(4-5)'
      do j = 1, this%wakes(i)%nsegments
        this%wakes(i)%segments(j)%connectivity  = arrayaux1(j, :)
        this%wakes(i)%segments(j)%adjacency     = arrayaux2(j, :)
        this%wakes(i)%segments(j)%indicesq(1:3) = (/ (r1, r1 = 3*(this%wakes(i)%segments(j)%connectivity(1)-1)+1,3*(this%wakes(i)%segments(j)%connectivity(1)-1)+3) /)
        this%wakes(i)%segments(j)%indicesq(4:6) = (/ (r1, r1 = 3*(this%wakes(i)%segments(j)%connectivity(2)-1)+1,3*(this%wakes(i)%segments(j)%connectivity(2)-1)+3) /)
        !print*, j, this%wakes(i)%segments(j)%connectivity, this%wakes(i)%segments(j)%adjacency
      end do

      !print*, '... Getting circulation rings: ring nr.(1); adj. seg.(2-5)'
      do j = 1, this%wakes(i)%nrings
        this%wakes(i)%rings(j)%adjacency = arrayaux3(j, :)
        !print*, j, this%wakes(i)%rings(j)%adjacency
      end do
      
      !< but only for the first row
      do j = 1, this%wakes(i)%nringss
        this%wakes(i)%rings(j)%indicesp = (/ (r1, r1 = 3*(j-1)+1,3*(j-1)+3) /)
        this%wakes(i)%rings(j)%indicesc = j
        
        j3a = 3*(this%wakes(i)%rings(j)%connectivity(1)-1)
        j3b = 3*(this%wakes(i)%rings(j)%connectivity(2)-1)
        j3c = 3*(this%wakes(i)%rings(j)%connectivity(3)-1)
        j3d = 3*(this%wakes(i)%rings(j)%connectivity(4)-1)

        !(/ (r1, r1 = 1,18) /)
        this%wakes(i)%rings(j)%indicesq(1:3)   = (/ (r1, r1 = j3a+1,j3a+3) /)
        this%wakes(i)%rings(j)%indicesq(4:6)   = (/ (r1, r1 = j3b+1,j3b+3) /)
        this%wakes(i)%rings(j)%indicesq(7:9)   = (/ (r1, r1 = j3c+1,j3c+3) /)
        this%wakes(i)%rings(j)%indicesq(10:12) = (/ (r1, r1 = j3d+1,j3d+3) /)
        
        !< indicesv and incicesq are related to the four nodes of ring: Mainly for determining control point ccordinates from nodal coordinates
        this%wakes(i)%rings(j)%indicesv(1:12)  = this%wakes(i)%rings(j)%indicesq(1:12)

        do k = 1,4
          this%wakes(i)%rings(j)%nodes(k)%node_aero%localid  = this%wakes(i)%rings(j)%connectivity(k)
          this%wakes(i)%rings(j)%nodes(k)%node_aero%globalid = j1ns + this%wakes(i)%rings(j)%connectivity(k)
        end do
      end do

      !< global indices for circulation at rings, coordinates and  velocity at nodes
      allocate(this%wakes(i)%indicesc(this%wakes(i)%nrings))
      this%wakes(i)%indicesc(1:this%wakes(i)%nrings) = j3cw + (/ (r1, r1 = 1,this%wakes(i)%nrings) /)

      allocate(this%wakes(i)%indicesp(3*this%wakes(i)%nrings))
      this%wakes(i)%indicesp(1:3*this%wakes(i)%nrings) = j3pw + (/ (r1, r1 = 1,3*this%wakes(i)%nrings) /)

      allocate(this%wakes(i)%indicesq(3*this%wakes(i)%nnodes))
      this%wakes(i)%indicesq(1:3*this%wakes(i)%nnodes) = j3qw + (/ (r1, r1 = 1,3*this%wakes(i)%nnodes) /)

      allocate(this%wakes(i)%indicesv(3*this%wakes(i)%nnodes))
      this%wakes(i)%indicesv(1:3*this%wakes(i)%nnodes) = j3qw + (/ (r1, r1 = 1,3*this%wakes(i)%nnodes) /)

      j3cw = j3cw + this%wakes(i)%nrings
      j3qw = j3qw + 3*this%wakes(i)%nnodes
      j3pw = j3pw + 3*this%wakes(i)%nrings
      j1nw = j1nw + this%wakes(i)%nnodes
      
      !< storing information for additional terms in linearization: pulling first wake row to bounded-vortex shet separation edge
      n_wake_edge_segments = 0
      if (allocated(this%wakes(i)%enodes)) n_wake_edge_segments = size(this%wakes(i)%enodes)-1
      do j = n_wake_edge_segments, 1, -1
        ID_first_row_wake_edge = n_wake_edge_segments-j+1
        call append_int_vec(temp_vec, this%wakes(i)%rings(ID_first_row_wake_edge)%adjacency(1), 1)
        call append_int_vec(temp_vec, this%wakes(i)%rings(ID_first_row_wake_edge)%adjacency(2), 1)
        call append_int_vec(temp_vec, this%wakes(i)%rings(ID_first_row_wake_edge)%adjacency(4), 1)
      end do
      
      !< marking first row wake segments (trailing edge), not to be considered in linearization of wakes segment velocity
      allocate(this%wakes(i)%segments_ID_first_row(2*n_wake_edge_segments+1))
      this%wakes(i)%segments_ID_first_row = temp_vec
      do j = 1, 2*n_wake_edge_segments+1
        ID_first_row_wake_edge = this%wakes(i)%segments_ID_first_row(j)
        this%wakes(i)%segments(ID_first_row_wake_edge)%first_row_wakeedge = 1
      end do
      
      !< identifiying first row (trailing edge) wake rings to be linearized
      allocate(this%wakes(i)%rings_ID_first_row_bounded_sheet(n_wake_edge_segments))
      this%wakes(i)%rings_ID_first_row = (/ (r1, r1 = 1, n_wake_edge_segments) /) 
      do j = n_wake_edge_segments, 1, -1
        this%wakes(i)%rings_ID_first_row_bounded_sheet(n_wake_edge_segments-j+1) = this%wakes(i)%erings(j)
      end do
      
      !< assign global bounded-vortex sheet nodesID for wake segment
      do j = 1, n_wake_edge_segments
        ! global bounded-vortex sheet nodesID for wake segment
        node1 = this%wakes(i)%enodes(n_wake_edge_segments-j+2)
        node2 = this%wakes(i)%enodes(n_wake_edge_segments-j+1)
        
        !< select segments in vortex ring
        segment1_number = this%wakes(i)%rings(j)%adjacency(1)
        segment2_number = this%wakes(i)%rings(j)%adjacency(2)
        segment3_number = this%wakes(i)%rings(j)%adjacency(3)
        segment4_number = this%wakes(i)%rings(j)%adjacency(4)
        
        this%wakes(i)%segments(segment1_number)%connectivity_to_bounded_sheet(1) = node1
        this%wakes(i)%segments(segment1_number)%connectivity_to_bounded_sheet(2) = node2
        
        this%wakes(i)%segments(segment2_number)%connectivity_to_bounded_sheet(1) = node2
        this%wakes(i)%segments(segment2_number)%connectivity_to_bounded_sheet(2) = 0
        
        this%wakes(i)%segments(segment3_number)%connectivity_to_bounded_sheet(1) = 0
        this%wakes(i)%segments(segment3_number)%connectivity_to_bounded_sheet(2) = 0
        
        this%wakes(i)%segments(segment4_number)%connectivity_to_bounded_sheet(1) = node1
        this%wakes(i)%segments(segment4_number)%connectivity_to_bounded_sheet(2) = 0
      end do
      
      deallocate(arrayaux1)
      deallocate(arrayaux2)
      deallocate(arrayaux3)

      !< determine transormation matrix for transformation between ring circulation and (global) segment circulation
      call this%wakes(i)%rings2segments_matrix()
    end do
  
!< allocating system vectors and coefficient Matrix A
    allocate(this%amatrix(this%tnrings, this%tnrings))
    allocate(this%inv_amatrix(this%tnrings, this%tnrings))
    allocate(this%rhs(this%tnrings))
    allocate(this%circulations(this%tnrings))
    allocate(this%oldcirculations(this%tnrings))
    allocate(this%deltaps(this%tnrings))
    allocate(this%fae(3*this%tnrings))
    allocate(this%faen(this%tncoordinates))
    allocate(this%areanormal(3*this%tnrings))
    
    allocate(this%fae_x(3*this%tnrings))
    allocate(this%fae_w(3*this%tnrings))
    allocate(this%fae_xn(this%tncoordinates))
    allocate(this%fae_wn(this%tncoordinates))
    
    allocate(this%deltav(3*this%tnrings))
    allocate(this%vmean(3*this%tnrings))
    allocate(this%vs(3*this%tnrings))
    allocate(this%vsv(3*this%tnrings))
    allocate(this%vwv(3*this%tnrings))

    !< nodal values as input
    allocate(this%qs_t(this%tncoordinates))
    allocate(this%vs_t(this%tnvelocities))

    !< nodal value of area normal, circulation and pressure, determine with interpolation
    allocate(this%circulations_nodal(this%tnnodes))
    allocate(this%deltaps_nodal(this%tnnodes))
    allocate(this%areanormal_nodal(this%tncoordinates))
    allocate(this%vmean_nodal(this%tnvelocities))

    !< arrays for linearization of non-penetration condition
    allocate(this%dAG (this%tnrings, 2*this%tncoordinates))
    allocate(this%drhs(this%tnrings, 2*this%tncoordinates))
    allocate(this%dG  (this%tnrings, 2*this%tncoordinates))
    allocate(this%dfae(this%tncoordinates, 2*this%tncoordinates))
    allocate(this%dfae_ring(3*this%tnrings, 2*this%tncoordinates))

    call sparse_initialize(this%sparse_dfae, this%tncoordinates, this%tncoordinates + this%tnvelocities)

    !< some solver buffer variables
    allocate(this%ipiv(this%tnrings))
    allocate(this%work(this%tnrings))

    this%amatrix(:, :)      = 0.0d0
    this%inv_amatrix(:, :)  = 0.0d0
    this%rhs(:)             = 0.0d0
    this%circulations(:)    = 0.0d0
    this%oldcirculations(:) = 0.0d0
    this%deltaps(:)         = 0.0d0
    this%fae(:)             = 0.0d0
    this%fae_x(:)           = 0.0d0
    this%fae_w(:)           = 0.0d0
    this%faen(:)            = 0.0d0
    this%fae_xn(:)          = 0.0d0
    this%fae_wn(:)          = 0.0d0
    this%areanormal(:)      = 0.0d0

    this%deltav(:)          = 0.0d0
    this%vmean(:)           = 0.0d0
    this%vs(:)              = 0.0d0
    this%vsv(:)             = 0.0d0
    this%vwv(:)             = 0.0d0

    this%vs_t(:)            = 0.0d0
    this%qs_t(:)            = 0.0d0

    this%circulations_nodal(:) = 0.0d0
    this%deltaps_nodal(:)      = 0.0d0
    this%areanormal_nodal(:)   = 0.0d0
    this%vmean_nodal(:)        = 0.0d0

    this%dAG  = 0.0d0
    this%drhs = 0.0d0
    this%dG   = 0.0d0
    this%dfae = 0.0d0
    this%dfae_ring = 0.0d0

    this%work = 0.0d0
    this%ipiv = 0

!< allocating linearized segment circulation array
    do i = 1, this%nsurfaces
      allocate(this%surfaces(i)%DGamma(this%surfaces(i)%nsegments,2*this%tncoordinates))
      allocate(this%surfaces(i)%DGamma_dense(this%surfaces(i)%nsegments,2*this%tncoordinates))
      allocate(this%surfaces(i)%DG(size(this%surfaces(i)%indicesc),2*this%tncoordinates))
    end do

!< initializing global vector of surfaces coordinates
    do i = 1, this%nsurfaces
        do j = 1, this%surfaces(i)%nnodes
          this%qs_t(this%surfaces(i)%indicesq(3*(j-1)+1:3*(j-1)+3)) = this%surfaces(i)%nodes(j)%q_0
        end do
    end do

    call this%actualize(.TRUE., .FALSE.)

  return
  end subroutine constructor

!< Subroutine for reading old simulation files and continuing
  subroutine read_oldsim(this)
    implicit none
    class(model_aero), intent(inout) :: this
    integer :: io_error, i, j, nsteps, tcounter, nvec1, nvec2, nvec3
    integer :: UnIn_t, UnIn_qs_nodal, UnIn_circs_cp, UnIn_wake_qw, UnIn_wake_vw, UnIn_wake_circw
    real(kind = 8), allocatable :: vec1(:), vec2(:), vec3(:), temp_circulations(:), temp_qs_t(:)
    real(kind = 8) :: time, real_tcounter
    character(len = 1024) :: char_t
    character(:), allocatable :: char_i

    logical :: bool_error = .FALSE., boolean_opened = .FALSE.

    ! initialize file identifier
    call GetNewUnit (UnIn_t)
    call GetNewUnit (UnIn_qs_nodal)
    call GetNewUnit (UnIn_circs_cp)
    call GetNewUnit (UnIn_wake_qw)
    call GetNewUnit (UnIn_wake_vw)
    call GetNewUnit (UnIn_wake_circw)

    !< reading result files from previous simulation
    call GetNewUnit (UnIn_t)
    open(unit = UnIn_t, file = this%oldsimfilename // '_uvlm_t.dres', status = 'old', iostat = io_error)
    if (io_error .ne. 0) bool_error = .TRUE.
    if (bool_error) go to 10

    !< determine number of rows in result files
    nsteps = 0
    read(UnIn_t, *, iostat = io_error)
    do
      if (io_error .ne. 0) exit
      read(UnIn_t, *, iostat = io_error)
      nsteps = nsteps + 1
    end do
    close(unit = UnIn_t)

    !< open result files
    call GetNewUnit (UnIn_t)
    open(unit = UnIn_t, file = this%oldsimfilename // '_uvlm_t.dres', status = 'old', iostat = io_error)
    if (io_error .ne. 0)  then
      bool_error = .TRUE.
      print*, '... cannot read ', this%oldsimfilename // '_uvlm_t.dres'
      go to 10
    end if

    call GetNewUnit (UnIn_qs_nodal)
    open(unit = UnIn_qs_nodal, file = this%oldsimfilename // '_uvlm_qs_nodal.dres', status = 'old', iostat = io_error)
    if (io_error .ne. 0)  then
      bool_error = .TRUE.
      print*, '... cannot read ', this%oldsimfilename // '_uvlm_qs_nodal.dres'
      go to 10
    end if

    call GetNewUnit (UnIn_circs_cp)
    open(unit = UnIn_circs_cp, file = this%oldsimfilename // '_uvlm_circs_cp.dres', status = 'old', iostat = io_error)
    if (io_error .ne. 0)  then
      bool_error = .TRUE.
      print*, '... cannot read ', this%oldsimfilename // '_uvlm_circs_cp.dres'
      go to 10
    end if

    !< jumping to the last last time step
    do i = 1, nsteps-1
      read(UnIn_t, *, iostat = io_error)
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error)  go to 10

      read(UnIn_qs_nodal, *, iostat = io_error)
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error) go to 10

      read(UnIn_circs_cp, *, iostat = io_error)
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error) go to 10

    end do

    !< reading results of the last time step
    read(UnIn_t, *, iostat = io_error) time, real_tcounter, this%delta_sim_time
    if (io_error .ne. 0) bool_error = .TRUE.
    if (bool_error) go to 10
    tcounter = IDNINT(real_tcounter)

    if (time .ge. this%totalt) then
      print*, 'error in simulation! Total simulation time lower equal then simulation time of previous steps!'
      stop
    end if

    allocate(temp_qs_t(size(this%qs_t)))
    allocate(temp_circulations(size(this%circulations)))

    read(UnIn_qs_nodal, *, iostat = io_error) temp_qs_t
    if (io_error .ne. 0) bool_error = .TRUE.
    if (bool_error) go to 10

    read(UnIn_circs_cp, *, iostat = io_error) temp_circulations
    if (io_error .ne. 0) bool_error = .TRUE.
    if (bool_error) go to 10

    !< assign surface results
    this%circulations(:)    = temp_circulations(:)
    this%oldcirculations(:) = temp_circulations(:)
    this%qs_t(:)            = temp_qs_t(:)

    this%time = time
    this%tcounter = tcounter + 1
    this%boolean_actualize_amatrix = .TRUE.
    this%nsteps = tcounter + IDNINT((this%totalt-this%time)/this%deltat) + 1

    do j = 1, this%nwakes
      write(char_t, '(i5)')  j
      allocate(char_i, source = trim(adjustl(char_t)))

      !< read position of wakes
      call GetNewUnit (UnIn_wake_qw)
      open(unit = UnIn_wake_qw, file = this%oldsimfilename // '_uvlm_qw' // char_i // '.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0)  then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_uvlm_qw' // char_i // '.dres'
        go to 10
      end if

      do i = 1, nsteps-1
        read(UnIn_wake_qw, *, iostat = io_error)
        if (io_error .ne. 0) bool_error = .TRUE.
        if (bool_error) go to 10
      end do

      call read_line_to_array(UnIn_wake_qw,vec1, 0, io_error)
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error) go to 10
      close(unit = UnIn_wake_qw)

      !< read velocity of wakes
      call GetNewUnit (UnIn_wake_vw)
      open(unit = UnIn_wake_vw, file = this%oldsimfilename // '_uvlm_vw' // char_i // '.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0)  then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_uvlm_vw' // char_i // '.dres'
        go to 10
      end if

      call read_line_to_array(UnIn_wake_vw,vec2, 0, io_error)
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error) go to 10
      close(unit = UnIn_wake_vw)

      !< read ring circulations of wakes
      call GetNewUnit (UnIn_wake_circw)
      open(unit = UnIn_wake_circw, file = this%oldsimfilename // '_uvlm_circw' // char_i // '.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0)  then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_uvlm_circw' // char_i // '.dres'
        go to 10
      end if

      call read_line_to_array(UnIn_wake_circw,vec3, 0, io_error)
      if (io_error .ne. 0) bool_error = .TRUE.
      if (bool_error) go to 10
      close(unit = UnIn_wake_circw)

      !< setting wake results
      nvec1 = size(vec1)
      nvec2 = size(vec2)
      nvec3 = size(vec3)

      if (nvec1 .eq. 0 .or. nvec2 .eq. 0 .or. nvec3 .eq. 0) bool_error = .TRUE.
      if (bool_error) go to 10

      if (nvec1 .ne. nvec2) bool_error = .TRUE.
      if (bool_error) go to 10

      !< determine number of wake rows to considered from result files and assigning wake positions and velocities
      if (nvec1 .le. this%wakes(j)%ncoordinates) then
        this%wakes(j)%q_t(1:size(vec1)) = vec1(:)
        this%wakes(j)%v_t(1:size(vec2)) = vec2(:)
        this%wakes(j)%tcounter = size(vec1)/3 / this%wakes(j)%nnodess - 1
      else
        this%wakes(j)%q_t(1:this%wakes(j)%ncoordinates) = vec1(1:this%wakes(j)%ncoordinates)
        this%wakes(j)%v_t(1:this%wakes(j)%ncoordinates) = vec2(1:this%wakes(j)%ncoordinates)
        this%wakes(j)%tcounter = this%wakes(j)%ncoordinates/3 / this%wakes(j)%nnodess - 1
      end if

      if (this%wakes(j)%tcounter .lt. tcounter+1) this%wakes(j)%bool_CutWake = .TRUE.
      if (this%wakes(j)%nrowstotal .lt. tcounter+1) this%wakes(j)%bool_CutWake = .TRUE.
      if (this%wakes(j)%tcounter .lt. this%wakes(j)%nstepstotal) this%wakes(j)%bool_CutWake = .FALSE.

      !< assigning wake circulation
      if (nvec3 .le. this%wakes(j)%nrings) then
        this%wakes(j)%ringscirculations(1:nvec3) = vec3(:)
      else
        this%wakes(j)%ringscirculations(1:this%wakes(j)%nrings) = vec3(1:this%wakes(j)%nrings)
      end if

      deallocate(char_i)
    end do

    print*, '... all aerodynamic simulation files read correctly'

    !< if error, then continue from here
10  continue

    inquire(unit=UnIn_qs_nodal, opened=boolean_opened)
    if (boolean_opened) close(unit = UnIn_qs_nodal)

    inquire(unit=UnIn_t, opened=boolean_opened)
    if (boolean_opened) close(unit = UnIn_t)

    inquire(unit=UnIn_circs_cp, opened=boolean_opened)
    if (boolean_opened) close(unit = UnIn_circs_cp)

    inquire(unit=UnIn_wake_qw, opened=boolean_opened)
    if (boolean_opened) close(unit = UnIn_wake_qw)

    inquire(unit=UnIn_wake_vw, opened=boolean_opened)
    if (boolean_opened) close(unit = UnIn_wake_vw)

    inquire(unit=UnIn_wake_circw, opened=boolean_opened)
    if (boolean_opened) close(unit = UnIn_wake_circw)

    if (bool_error) then
      print*, 'warning in aero simulation! Error in files! Simulation restarts!'
      this%boolean_oldsim = .FALSE.
    end if

    return
  end subroutine read_oldsim

!< Subroutine to read wake from file
  subroutine read_oldwake(this)
    implicit none
    class(model_aero), intent(inout) :: this
    integer :: io_error, i, j, nsteps, tcounter, nvec1, nvec2, nvec3, UnIn
    real(kind = 8), allocatable :: vec1(:), vec2(:), vec3(:), temp_circulations(:), temp_qs_t(:)
    real(kind = 8) :: time, real_tcounter
    character(len = 1024) :: char_t
    character(:), allocatable :: char_i

    logical :: bool_error = .FALSE.

    !< reading result files from previous simulation
    do j = 1, this%nwakes
      nsteps = 0
      write(char_t, '(i5)')  j
      allocate(char_i, source = trim(adjustl(char_t)))

      call GetNewUnit (UnIn)
      open(unit = UnIn, file = this%oldsimfilename // '_uvlm_qw' // char_i // '.dres', status = 'old', iostat = io_error)
      if (io_error .eq. 0) then
        do ! loop over all row until end of file
          read(UnIn, *, iostat = io_error)
          if (io_error .lt. 0) exit
          nsteps = nsteps + 1
        end do
        close(unit = UnIn)

        !< read position of wakes
        open(unit = UnIn, file = this%oldsimfilename // '_uvlm_qw' // char_i // '.dres', status = 'old', iostat = io_error)
        do i = 1, nsteps - 1
          read(UnIn, *, iostat = io_error)
        end do
        call read_line_to_array(UnIn,vec1, 0, io_error)
        close(unit = UnIn)
        if (io_error .ne. 0) bool_error = .TRUE.
      else
        print*,'... cannot read ', this%oldsimfilename // '_uvlm_qw' // char_i // '.dres'
        bool_error = .TRUE.
      end if

      call GetNewUnit (UnIn)
      open(unit = UnIn, file = this%oldsimfilename // '_uvlm_circw' // char_i // '.dres', status = 'old', iostat = io_error)
      if (io_error .ne. 0)  then
        bool_error = .TRUE.
        print*, '... cannot read ', this%oldsimfilename // '_uvlm_circw' // char_i // '.dres'
      else
        !< read ring circulations of wakes
        call read_line_to_array(UnIn,vec3, 0, io_error)
        close(unit = UnIn)
        if (io_error .ne. 0) bool_error = .TRUE.
      end if

      deallocate(char_i)
      if (bool_error) go to 10

      !< setting wake results
      nvec1 = size(vec1)
      nvec3 = size(vec3)

      this%wakes(j)%q_t(1:nvec1) = vec1(:)
      this%wakes(j)%v_t(1:nvec1) = 0.0d0
      this%wakes(j)%ringscirculations(1:nvec3) = vec3(:)
      this%wakes(j)%tcounter = nvec1/3 / this%wakes(j)%nnodess - 1
      this%nsteps = nsteps

    end do

    !< if error, then continue from here
10    continue

    if (bool_error) then
      print*, 'Error in reading wake files! No result or error in files found! Simulation aborts!'
      this%boolean_abort = .TRUE.
    end if

    return
  end subroutine read_oldwake

!< Subroutine aero solver
  subroutine solver(this)
    implicit none
    class(model_aero), intent(inout) :: this
    integer :: info
    real(kind = 8), allocatable :: amatrix(:,:)

    integer :: i, j

    allocate(amatrix(this%tnrings,this%tnrings))
    amatrix = 0.0d0

    !< assign surface geometry of bounded vortex sheets
    if (allocated(this%arrqs_t)) this%qs_t(:) = this%arrqs_t(this%tcounter+1,1:this%tncoordinates)
    if (allocated(this%arrvs_t)) this%vs_t(:) = this%arrvs_t(this%tcounter+1,1:this%tnvelocities)

    !< set initial free-field velocity for external flow field; for prescirbed flow-field from file, this is valid for outside of flow-field grid
    this%airflow%vinfinity = this%airflow%dir * this%evaluateamplitude(this%time)

    !< update geometry and velocity bounded vortex sheets
    call this%actualize(.TRUE., .FALSE.)
    this%oldcirculations = this%circulations

    !< update the wake geometry
    if (this%tcounter  .eq. 0) then
      call this%wake_init()
    else
      call this%update_wake()
    end if

    !< A, rhs and solve for g
    call this%evaluaterhs(this%cutoff)
    if (this%bool_const_a_matrix) then
      if (this%boolean_actualize_amatrix) then
        call this%evaluateamatrix(this%cutoff)
        amatrix = this%amatrix
        call dgetrf(this%tnrings, this%tnrings, amatrix, this%tnrings, this%ipiv, info )             !< LU decomposition and inverse of A
        call dgetri(this%tnrings, amatrix, this%tnrings, this%ipiv, this%work, this%tnrings, info )  !< inverse of A
        this%inv_amatrix = amatrix
        this%boolean_actualize_amatrix = .FALSE.
      end if
      this%circulations = matmul(this%inv_amatrix,this%rhs)
    else
      call this%evaluateamatrix(this%cutoff)
      amatrix = this%amatrix
      this%circulations = this%rhs
      call dgesv(this%tnrings, 1, amatrix, this%tnrings, this%ipiv, this%circulations, this%tnrings, info)
    end if

    !< actualize ring intensity, evaluate loads, calculate wake veloctiy
    call this%actualize(.FALSE., .TRUE.)
    call this%modelfae()

    !< calculate sensitivity matrices for last simulation step
    !if (this%tcounter+1 .eq. this%nsteps) then
    !  call this%modeldfae(.True.)
    !end if

    !< update velocities and old circulation in wake nodes
    call this%wake_velocity(this%cutoff)

  return
  end subroutine solver

!< solver for one step calculation of aerodynamic forces considering prescribed (from file) wake
subroutine solver_num_sensitivity(this)
  implicit none

  class(model_aero), intent(inout) :: this

  integer :: info, i, io_error, nsteps, UnIn
  real(kind = 8), allocatable :: amatrix(:,:)

  !< read file for old circulations
  call GetNewUnit (UnIn)
  open(unit = UnIn, file = this%oldsimfilename // '_uvlm_circs_cp.dres', status = 'old', iostat = io_error)
  if (io_error .ne. 0)  then
    print*, '... cannot read ', this%oldsimfilename // '_uvlm_circs_cp.dres'
    stop
  end if

  !< jumping to the before last time step
  do i = 1, this%nsteps-2
    read(UnIn, *, iostat = io_error)
    if (io_error .ne. 0) exit
  end do

  read(UnIn, *, iostat = io_error) this%oldcirculations(:)
  if (io_error .ne. 0) stop
  close(unit = UnIn)

  this%tcounter = this%nsteps - 1

  allocate(amatrix(this%tnrings,this%tnrings))
  amatrix = 0.0d0

  !< update wake geometry and properties
  do i = 1,this%nwakes
    call this%wakes(i)%setproperty()
  end do

  !< set initial free-field velocity for external flow field; for prescirbed flow-field from file, this is valid for outside of flow-field grid
  this%airflow%vinfinity = this%airflow%dir * this%evaluateamplitude(this%time)

  !< assign surface geometry of bounded vortex sheets
  if (allocated(this%arrqs_t)) this%qs_t(:) = this%arrqs_t(this%tcounter+1,1:this%tncoordinates)

  !< update geometry and velocity bounded vortex sheets
  call this%actualize(.TRUE., .FALSE.)

  !< A, rhs and solve for g
  call this%evaluaterhs(this%cutoff)
  call this%evaluateamatrix(this%cutoff)
  amatrix = this%amatrix
  call dgetrf(this%tnrings, this%tnrings, amatrix, this%tnrings, this%ipiv, info )             !< LU decomposition and inverse of A
  call dgetri(this%tnrings, amatrix, this%tnrings, this%ipiv, this%work, this%tnrings, info )  !< inverse of A
  this%inv_amatrix = amatrix
  this%circulations = matmul(this%inv_amatrix,this%rhs)

  !< actualize ring intensity, evaluate loads, calculate wake veloctiy
  call this%actualize(.FALSE., .TRUE.)
  call this%modelfae()
  call this%modeldfae(.True.)

  deallocate(amatrix)
  return
end subroutine solver_num_sensitivity

!< Subroutine to read input files
  subroutine readinginputs(this)
    use my_constants_aero, only : eps_val
    implicit none

    class(model_aero), intent(inout) :: this

    integer :: io_error_surfaces, io_error_wakes, io_error_simulation, io_error
    integer :: i, j
    integer, parameter :: linejump = 3
    integer, allocatable :: edge(:, :), wake_prop(:,:)
    real(kind = 8), allocatable :: wake_diff(:)
    real(kind = 8) :: time
    character(100) :: vinf_sort, stroldfilename, strfilename, strType, strfilename_windfield

!< reading input file of vortex sheet
    print*, 'Parsing bounded vortex elements'
    open(unit = 5000, file = 'surfaceinput.txt', status = 'old', iostat = io_error_surfaces)
    if (io_error_surfaces == 0) then !< if the file surfaceinput exists, then read the information.
       do i = 1, linejump
          read(5000, *)
       end do
       read(5000, *) this%nsurfaces
       allocate(this%surfaces(this%nsurfaces))
       do i = 1, this%nsurfaces

          do j = 1, linejump
             read(5000, *)
          end do

          read(5000, *) this%surfaces(i)%nnodes, this%surfaces(i)%nrings, this%surfaces(i)%nnodes1, this%surfaces(i)%nnodes2
          this%surfaces(i)%ncoordinates = 3*this%surfaces(i)%nnodes

          allocate(this%surfaces(i)%nodes(this%surfaces(i)%nnodes))
          allocate(this%surfaces(i)%rings(this%surfaces(i)%nrings))

          do j = 1, linejump
             read(5000, *)
          end do

          do j = 1, this%surfaces(i)%nnodes
             read(5000, *) this%surfaces(i)%nodes(j)%q_0(:)
          end do

          do j = 1, linejump
             read(5000, *)
          end do

          do j = 1, this%surfaces(i)%nrings
             read(5000, *) this%surfaces(i)%rings(j)%connectivity
          end do
       end do
    else
      this%boolean_abort = .TRUE.
      print*, 'Warning: no surfaceinput.txt found'
    end if
    close(unit = 5000)

!< reading wake information
    print*, 'Parsing unbounded vortex elements'
    open(unit = 5001, file = 'wakeinput.txt', status = 'old', iostat = io_error_wakes)
    if (io_error_wakes == 0) then ! if the file wakeinput exists, then read the information.

      do i = 1, linejump
         read(5001, *)
      end do

      read(5001, *) this%nwakes, this%nwakesproperties
      allocate(this%wakes(this%nwakes))
      do i = 1, this%nwakes
        do j = 1, linejump
          read(5001, *)
        end do

        read(5001, *) this%wakes(i)%surface, this%wakes(i)%nringss, this%wakes(i)%nproperty
        allocate(edge(this%wakes(i)%nringss, 3))
        do j = 1, linejump
          read(5001, *)
        end do

        !< Getting wake lines and corresponding circulation ring
        do j = 1, this%wakes(i)%nringss
          read(5001, *) edge(j, 1), edge(j, 2), edge(j, 3)
        end do
        this%wakes(i)%nnodess = this%wakes(i)%nringss+1

        allocate(this%wakes(i)%erings(this%wakes(i)%nringss))
        allocate(this%wakes(i)%enodes(this%wakes(i)%nnodess))

        this%wakes(i)%erings(1:this%wakes(i)%nringss) = edge(:, 3)
        this%wakes(i)%enodes(1:this%wakes(i)%nringss) = edge(:, 1)
        this%wakes(i)%enodes(  this%wakes(i)%nnodess) = edge(this%wakes(i)%nringss, 2)
        deallocate(edge)
      end do

      do j = 1, linejump
        read(5001, *)
      end do

      allocate(wake_prop(this%nwakesproperties,2))
      wake_prop(:,:) = 0
      allocate(wake_diff(this%nwakesproperties))
      wake_diff(:) = 0.0d0
      do i = 1, this%nwakesproperties
        read(5001, *, iostat = io_error) wake_prop(i,1), wake_prop(i,2), wake_diff(i)
      end do

      !< assigning wake properties to wakes
      do i = 1, this%nwakes
        this%wakes(i)%nrowstotal            = wake_prop(this%wakes(i)%nproperty,1)
        this%wakes(i)%nrowstoconsider       = wake_prop(this%wakes(i)%nproperty,2)
        this%wakes(i)%diffusion_coefficient = wake_diff(this%wakes(i)%nproperty)
      end do
    else
      print*, 'Warning: no wakeinput.txt found'
    end if
    close(unit = 5001)

!< reading simulation information
    print*, 'Parsing aero-simulation input file'
    open(unit = 5002, file = 'simulationinput_aero.txt', status = 'old', iostat = io_error_simulation)
    if (io_error_simulation == 0) then ! if the file simulationinput exists, then read the information.

      do i = 1, linejump
        read(5002, *)
      end do
      !< reading filename
      read(5002, *) strfilename, stroldfilename
      allocate(this%output_filename, source =  trim(adjustl(strfilename)))
      if (trim(adjustl(stroldfilename)) .ne. 'none' .and. trim(adjustl(stroldfilename)) .ne. '!!') then
        this%boolean_oldsim = .TRUE.
        allocate(this%oldsimfilename, source =  trim(adjustl(stroldfilename)))
      end if
      if (trim(adjustl(stroldfilename)) .eq. '!!') backspace(5002)

      !< reading simulation settings
      do i = 1, linejump
          read(5002, *)
      end do
      read(5002, *) this%totalt, this%deltat, this%cutoff, strType
      if (trim(adjustl(strType)) .eq. '!!') then
        backspace(5002)
        strType = 'aero'
      end if
      allocate(this%strSimuType, source =  trim(adjustl(strType)))

      !< determine number of steps - due to numerical reasons of double precision in this way
      this%nsteps = 0
      time = 0.0d0
      do while (time .le. this%totalt)
        time = time + real(this%deltat)
        this%nsteps = this%nsteps + 1
      end do
      this%nsteps = IDNINT(this%totalt/this%deltat) + 1

      do i = 1, linejump
        read(5002, *)
      end do

      !< reading wind settings
      read(5002, *) vinf_sort, this%airflow%density, this%airflow%intensity, this%airflow%duration, this%airflow%dir
      allocate(this%airflow%sort, source =  trim(adjustl(vinf_sort)))

      !< if sort of airflow is file then read airflow informations
      if(adjustl(trim(this%airflow%sort)) .eq. 'file') then
        do i = 1, linejump
          read(5002, *)
        end do

        read(5002, *) strfilename_windfield
        allocate(this%airflow%strfilename_windfield, source = adjustl(trim(strfilename_windfield)))

        do i = 1, linejump
          read(5002, *)
        end do
        read(5002, *) this%airflow%center
      else
        !< checking for correct airflow sort - if not  correct defined, then boolean_abort_airflow = .TRUE.
        this%airflow%vinfinity = this%airflow%dir * this%evaluateamplitude(this%time)
        if (this%boolean_abort_airflow) then
          print*,'Warning: continuing with constant flow field...'
          this%airflow%sort = 'constant'
          this%boolean_abort = .FALSE.
        end if
      end if

    else
      this%boolean_abort = .TRUE.
      print*,'Warning: no simulation file for DeSiO-Aero available...'
    end if
    close(unit = 5002)

    !< read airflow informations, if available
    if (allocated(this%airflow%strfilename_windfield)) then
      this%airflow%sort = 'constant'
      call this%airflow%getairflowDatafromfile()
      if (this%airflow%boolean_abort) then
        this%boolean_abort = .TRUE.
        return
      end if
    end if

    this%boolean_model_read = .TRUE.
  return
  end subroutine readinginputs

!< Subroutine to set and actualize vortex sheeet properties (geometry and circulation)
  subroutine actualize(this, boolean_q, boolean_circ)
    implicit none

    class(model_aero), intent(inout) :: this
    logical, intent(in) :: boolean_q, boolean_circ
    integer :: i

    do i = 1, this%nsurfaces
      if (boolean_circ) then
        this%surfaces(i)%ringscirculations    = this%circulations(this%surfaces(i)%indicesc)
        this%surfaces(i)%oldringscirculations = this%oldcirculations(this%surfaces(i)%indicesc)
        call this%surfaces(i)%setproperty()
      end if

      if (boolean_q) then
        this%surfaces(i)%q_t = this%qs_t(this%surfaces(i)%indicesq)
        this%surfaces(i)%v_t = this%vs_t(this%surfaces(i)%indicesv)
        call this%surfaces(i)%setkinematic()
      end if
    end do

    return
  end subroutine actualize

!< Subrouinte to evaluate A-Matrix, which depends on the radius between targetpoint (on which velocity is calculated) and circulation (segment) points.
  subroutine evaluateamatrix(this, cutoff_opt)
    implicit none

    class(model_aero), intent(inout) :: this
    integer i, j, k, l, indexk, indexl
    real(kind = 8) :: circulation_temp
    real(kind = 8), optional :: cutoff_opt

    indexk = 0
    do i = 1, this%nsurfaces
       do k = 1, this%surfaces(i)%nrings
          indexl = 0
          do j = 1, this%nsurfaces
!!$OMP PARALLEL DO &
!!$OMP SHARED(i, k, j, indexk, indexl, cutoff_opt) &
!!$OMP PRIVATE(l, circulation_temp)
            do l = 1, this%surfaces(j)%nrings
              !< setting unit ring circulation for A
              circulation_temp = this%surfaces(j)%rings(l)%circulation
              this%surfaces(j)%rings(l)%circulation = 1.0d0
              !< setting ring cutoff
              this%surfaces(j)%rings(l)%cutoff = 0.0d0
              if (present(cutoff_opt)) this%surfaces(j)%rings(l)%cutoff = cutoff_opt
              !< storing velocitiy due to unit ring circulations for post processing, i.e., pressure difference
              this%surfaces(i)%rings(k)%avector(3*(indexl+l-1)+1:3*(indexl+l-1)+3) = this%surfaces(j)%rings(l)%velocity(this%surfaces(i)%rings(k)%center)
              !< calculating coeeficient matrix A
              this%amatrix(indexk+k, indexl+l) = dot_product(this%surfaces(i)%rings(k)%avector(3*(indexl+l-1)+1:3*(indexl+l-1)+3),this%surfaces(i)%rings(k)%normal)
              !< assigning old circulation to rings circulation
              this%surfaces(j)%rings(l)%circulation = circulation_temp
            end do
!!$OMP  END PARALLEL DO
            indexl = indexl+this%surfaces(j)%nrings
          end do
       end do
       indexk = indexk+this%surfaces(i)%nrings
    end do

    return
  end subroutine evaluateamatrix

!< subroutine to calculate right-hand side of system of equations: wake induced velocity + vinf
  subroutine evaluaterhs(this, cutoff_opt)
    implicit none

    class(model_aero), intent(inout) :: this
    integer i, j, indexj, k
    real(kind = 8) :: vwakes(3), cutoff
    real(kind = 8), optional :: cutoff_opt

    cutoff = 0.0d0
    if (present(cutoff_opt)) cutoff = cutoff_opt
    
    indexj = 0
    do i = 1, this%nsurfaces
      !$OMP PARALLEL DO &
      !$OMP SHARED(i, indexj, cutoff) &
      !$OMP PRIVATE(j, k, vwakes)
       do j = 1, this%surfaces(i)%nrings
          vwakes(:) = 0.0d0
          do k = 1, this%nwakes
            if (allocated(this%wakes(k)%segments_t)) then
              vwakes = vwakes+this%wakes(k)%velocity(this%surfaces(i)%rings(j)%center, cutoff, this%time)
            end if
          end do

          !< storing vwake at ring for load calculation in evaluateloads
          this%surfaces(i)%rings(j)%vwakes = vwakes

          !< update free field velocity for each ring
          if (this%airflow%boolean_actualize) then
            if (this%airflow%boolean_external) then
              this%surfaces(i)%rings(j)%vinfinity = this%airflow%airflowvinfinity(this%surfaces(i)%rings(j)%center, this%time)
            else
              this%surfaces(i)%rings(j)%vinfinity = this%airflow%vinfinity
            end if
          end if

          !< right-hand side of system of equations
          this%rhs(indexj+j) =-dot_product(this%surfaces(i)%rings(j)%vinfinity + this%surfaces(i)%rings(j)%vwakes - this%surfaces(i)%rings(j)%vsurf, this%surfaces(i)%rings(j)%normal)
       end do
      !$OMP  END PARALLEL DO
       indexj = indexj+this%surfaces(i)%nrings
    end do
    return
  end subroutine evaluaterhs

!< Subroutine to update first wake row
  subroutine wake_update_first_row(this, boolean_geometry_opt, boolean_circulation_opt)
    implicit none

    class(model_aero), intent(inout) :: this
    integer :: i, j
    logical :: boolean_geometry = .True. , boolean_circulation = .True.
    logical, optional :: boolean_geometry_opt, boolean_circulation_opt
    
    boolean_geometry = .True.
    if (present(boolean_geometry_opt)) boolean_geometry = boolean_geometry_opt
    
    boolean_circulation = .True.
    if (present(boolean_circulation_opt)) boolean_circulation = boolean_circulation_opt
      
    !< set first row of wake
    do i = 1, this%nwakes
          !< update wake property should come before setgeometry, as it defines new wake vortex segment
          if (boolean_circulation) then
            if (this%tcounter .gt. 0) then
              this%wakes(i)%ringscirculations(1:this%wakes(i)%nringss) = this%surfaces(this%wakes(i)%surface)%rings(this%wakes(i)%rings_ID_first_row_bounded_sheet)%circulation
              call this%wakes(i)%setproperty()
            end if
          end if
        !< update wake property should come first, as it defines new wake vortex segment
          if (boolean_geometry) then
            do j = 1, this%wakes(i)%nnodess
              this%wakes(i)%q_t(this%wakes(i)%nodes(j)%indicesq) = this%surfaces(this%wakes(i)%surface)%nodes(this%wakes(i)%enodes(this%wakes(i)%nnodess-(j-1)))%q_t
            end do
            call this%wakes(i)%setgeometry(.True.)
          end if
     end do

    return
  end subroutine wake_update_first_row

!< Subroutine to initialize the wake first row
  subroutine wake_init(this)
    implicit none

    class(model_aero), intent(inout) :: this
    integer :: i, j

  !< set wake first row as initialization
    do i = 1, this%nwakes
        do j = 1, this%wakes(i)%nnodess
          this%wakes(i)%q_t(this%wakes(i)%nodes(j)%indicesq) = this%surfaces(this%wakes(i)%surface)%nodes(this%wakes(i)%enodes(this%wakes(i)%nnodess-(j-1)))%q_t
        end do
    end do

    return
  end subroutine wake_init

!< Subroutine to calculate velocities in wake nodes
 subroutine wake_velocity(this, cutoff_opt)
    implicit none
    class(model_aero), intent(inout) :: this
    integer :: i, j, l
    real(kind = 8) :: vsurfaces(3), vwakes(3), cutoff
    real(kind = 8), optional :: cutoff_opt
    real(kind = 8) :: temp_q_t_wakes(3)

    cutoff = 0.0d0
    if (present(cutoff_opt)) cutoff = cutoff_opt

    !< loop over all wakes
    do i = 1, this%nwakes

      !$OMP PARALLEL DO &
      !$OMP SHARED(i) &
      !$OMP PRIVATE(l, j, vwakes, vsurfaces)
      !< loop over all wake nodes
      do j = 1, this%wakes(i)%nnodess*(this%wakes(i)%tcounter+1)

          !< influence of surface cirulations to wake velocity at wake nodes
          vsurfaces(:) = 0.0d0
          do l = 1, this%nsurfaces
            temp_q_t_wakes=this%wakes(i)%q_t(this%wakes(i)%nodes(j)%indicesq)
            vsurfaces = vsurfaces+this%surfaces(l)%velocity(temp_q_t_wakes, cutoff)
            !vsurfaces = vsurfaces+this%surfaces(l)%velocity(this%wakes(i)%q_t(this%wakes(i)%nodes(j)%indicesq), cutoff)
          end do

          !< influence of wake cirulations to wake velocity at wake node
          vwakes(:) = 0.0d0
          if (this%tcounter > 0) then
             do l = 1, this%nwakes
                temp_q_t_wakes=this%wakes(i)%q_t(this%wakes(i)%nodes(j)%indicesq)
                !vwakes = vwakes+this%wakes(l)%velocity(this%wakes(i)%q_t(this%wakes(i)%nodes(j)%indicesq), cutoff, this%time)
                vwakes = vwakes+this%wakes(l)%velocity(temp_q_t_wakes, cutoff, this%time)
             end do
          end if

          !< update free-field velocity
          if (this%airflow%boolean_actualize) then
            if (this%airflow%boolean_external) then
              this%wakes(i)%nodes(j)%vinfinity = this%airflow%airflowvinfinity(this%wakes(i)%q_t(this%wakes(i)%nodes(j)%indicesq), this%time)
            else
              this%wakes(i)%nodes(j)%vinfinity = this%airflow%vinfinity
            end if
          end if

          !< calculating veloctiy in wake nodes
          this%wakes(i)%v_t(this%wakes(i)%nodes(j)%indicesq) = this%wakes(i)%nodes(j)%vinfinity + vwakes + vsurfaces
      end do
      !$OMP  END PARALLEL DO
    end do

    return
  end subroutine wake_velocity

!< Subroutine to calculate convection of wake nodes
  subroutine update_wake(this)
    implicit none

    class(model_aero), intent(inout) :: this
    real(kind = 8), allocatable :: temp_ringscirculations(:)
    integer :: i, j, k3
    integer :: nnodes, nnodes_old

!< Actualization of nodal position and circulations for the current load step
    do i = 1, this%nwakes

      !< update time counter for wake
      if (this%wakes(i)%bool_CutWake .eqv. .FALSE.) then
        if (this%wakes(i)%tcounter .lt. this%wakes(i)%nstepstotal) this%wakes(i)%tcounter = this%wakes(i)%tcounter + 1 ! this%tcounter
      else
        this%wakes(i)%tcounter = this%wakes(i)%nrowstoconsider
        !< reallocate arrays: only one time to do
        if (this%wakes(i)%bool_reallocate .eqv. .TRUE.) then
          call this%wakes(i)%reallocate_arrays()
        end if
      end if

! Actualization of the nodal positions
      nnodes       = this%wakes(i)%nnodess*(this%wakes(i)%tcounter+1)
      nnodes_old   = this%wakes(i)%nnodess*(this%wakes(i)%tcounter)
      this%wakes(i)%q_t(3*this%wakes(i)%nnodess+1:3*nnodes) = this%wakes(i)%q_t(1:3*nnodes_old) + this%wakes(i)%v_t(1:3*nnodes_old)*this%deltat

! Actualization of the first wake row, according to updated surface kinematic
      do j = 1, this%wakes(i)%nnodess
        this%wakes(i)%q_t_first_row(this%wakes(i)%nodes(j)%indicesq) = this%surfaces(this%wakes(i)%surface)%nodes(this%wakes(i)%enodes(this%wakes(i)%nnodess-(j-1)))%q_t
      end do
      this%wakes(i)%q_t(1:3*this%wakes(i)%nnodess) = this%wakes(i)%q_t_first_row(:)

! Actualization of the circulation
      this%wakes(i)%temp_ringscirculations(1:this%wakes(i)%nringss*(this%wakes(i)%tcounter-1)) = this%wakes(i)%ringscirculations(1:this%wakes(i)%nringss*(this%wakes(i)%tcounter-1))
      !$OMP PARALLEL DO &
      !$OMP SHARED(i) &
      !$OMP PRIVATE(j)
      do j = this%wakes(i)%nringss*(this%wakes(i)%tcounter-1), 1, -1
        this%wakes(i)%ringscirculations(j+this%wakes(i)%nringss) =  this%wakes(i)%temp_ringscirculations(j)
      end do
      !$OMP END PARALLEL DO

      !< assigning surface circulations at the wake edge to first row of wake circulations
      do j = 1, this%wakes(i)%nringss
         this%wakes(i)%ringscirculations(this%wakes(i)%nringss+1-j) = this%surfaces(this%wakes(i)%surface)%rings(this%wakes(i)%erings(j))%circulation
      end do

      !< updating segment and ring circulations and segment and ring geometry.
      !< important is that setproperty() should come first, as in this routine the new segments are defined
      call this%wakes(i)%setproperty()
      call this%wakes(i)%setgeometry()

      !< set flag for cutting the wake and reallocate the arrays of the wakes
      if (this%wakes(i)%tcounter .eq. this%wakes(i)%nstepstotal) then
        this%wakes(i)%bool_CutWake = .TRUE.
      end if
      if (this%wakes(i)%bool_CutWake .eqv. .TRUE.) this%wakes(i)%tcounter = this%wakes(i)%nrowstoconsider

    end do

    return
  end subroutine update_wake

! ==================================================================================================================================================
!< Subroutine to calculate aerodynamic force vector and tangential matrix
  subroutine modeldfae(this, boolean_kae)
    implicit none

    class(model_aero), intent(inout) :: this
    integer :: i, j, k, l, m, i1
    integer :: indicessegmentcirc(4)
    integer :: indicesq_rings_ij(12), indicesv_rings_ij(12), indicesq_rings_kl(12), indicesq_nodes_kl(3), indicesv_nodes_kl(3)
    integer :: indicesp_rings_ij(3)
    real(kind = 8) :: kg1(3,12), factor_tvec(12)
    real(kind = 8) :: kvs1(3,3), kvs2(3,12), kdv1(3,12), kvsv_G_q(3,3), kvsv_G_v(3,3), kdv_dcirc_q(3,3), kdv_dcirc_v(3,3)
    real(kind = 8), dimension(3,12) :: kg_a, kg_n, kg_Dv_lvec, kg_vwakes, kg_vsurf
    real(kind = 8), dimension(3,3) :: k33_vsv_DG, k33_Dv_DG, k33_G_DG, kg_vwakes_seg(3,3)
    real(kind = 8) :: knn_temp(12,12), k123_temp(12,3)
    logical, intent(in) :: boolean_kae
    real(kind = 8), allocatable :: dfae_ring(:,:)
    logical :: boolean_output_dfae_ring = .FALSE.
    type(tangent_matrices3_3x3):: kvwv_matrices
    integer :: indicesq_nodes_ij_l_1(3), indicesq_nodes_ij_l_2(3)
    integer :: node1, node2, segment_ID
    
    this%dfae = 0.0d0

    if (boolean_output_dfae_ring) then
      allocate(dfae_ring(3*this%tnrings,2*this%tncoordinates))
      dfae_ring = 0.0d0
    end if

    kg_a = 0.0d0
    kg_n = 0.0d0
    kg_Dv_lvec = 0.0d0
    kg_vwakes = 0.0d0
    kg_vsurf = 0.0d0
    k33_vsv_DG = 0.0d0
    k33_Dv_DG = 0.0d0
    k33_G_DG = 0.0d0
    kg_vwakes_seg = 0.0d0

    ! computing tangent matrices of aerodynamical forces - linearization of Bernoulli's equation
    call this%linearized_non_penetration_condition()

    if (boolean_kae) then
      i1 = 0
      do i = 1, this%nsurfaces

          do j = 1, this%surfaces(i)%nrings

            if (this%surfaces(i)%boolean_lifting_surface .eqv. .FALSE.) then
              exit
            end if

            indicesq_rings_ij = this%surfaces(i)%indicesq(this%surfaces(i)%rings(j)%indicesq)
            indicesv_rings_ij = this%tncoordinates + indicesq_rings_ij
            indicesp_rings_ij = this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp)

            !< linearization of surface geometry (normal and area)
            kg_a = this%surfaces(i)%rings(j)%deltaps * outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%karea)
            kg_n = this%surfaces(i)%rings(j)%deltaps * this%surfaces(i)%rings(j)%area * this%surfaces(i)%rings(j)%knormal

            !< linearization of Dv: derivative with respect to ring nodes of ring j, linearization with constant circulation
            call this%surfaces(i)%rings(j)%tangent_deltav(this%surfaces(i)%segmentscirculations, this%surfaces(i)%oldsegmentscirculations, kdv1, indicessegmentcirc, factor_tvec)
            kg_Dv_lvec = this%airflow%density * this%surfaces(i)%rings(j)%area * &
                            outer(this%surfaces(i)%rings(j)%normal,matmul(this%surfaces(i)%rings(j)%vmean-this%surfaces(i)%rings(j)%vsurf,kdv1))

            !< linearization of vwakes with respect to bounded-vortex sheet i ring j: kvwakes is calculated in the routine linearized_non_penetration_condition()
            kg_vwakes = this%airflow%density * this%surfaces(i)%rings(j)%area * matmul(outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%deltav), this%surfaces(i)%rings(j)%kvwakes)

            !< linearization of segments coordinates with repsect to nodal coordinates on separation edge; additional terms due to updated first wake row
            if (this%tcounter .ne. 0 .and. this%boolean_lin_wake_row .eqv. .TRUE.) then
              do k = 1, this%nwakes
                do l = 1, size(this%wakes(k)%segments_ID_first_row)
                  segment_ID = this%wakes(k)%segments_ID_first_row(l)
                  kvwv_matrices = this%wakes(k)%segments(segment_ID)%tangent_velocity(this%surfaces(i)%rings(j)%center)
              
                  !< assining to global nodes of bounded-vortex sheet
                  node1 = this%wakes(k)%segments(segment_ID)%connectivity_to_bounded_sheet(1)
                  node2 = this%wakes(k)%segments(segment_ID)%connectivity_to_bounded_sheet(2)
                  if (node1 .ne. 0) then
                    indicesq_nodes_ij_l_1 = this%surfaces(this%wakes(k)%surface)%indicesq(this%surfaces(this%wakes(k)%surface)%nodes(node1)%indicesq)
            
                    kg_vwakes_seg = this%airflow%density * this%surfaces(i)%rings(j)%area * &
                        matmul(outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%deltav), kvwv_matrices%kvs2)                  
                  
                    k123_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),kg_vwakes_seg)
                    this%dfae(indicesq_rings_ij,indicesq_nodes_ij_l_1) = this%dfae(indicesq_rings_ij,indicesq_nodes_ij_l_1) + k123_temp
                    if (boolean_output_dfae_ring) dfae_ring(indicesp_rings_ij,indicesq_nodes_ij_l_1) = dfae_ring(indicesp_rings_ij,indicesq_nodes_ij_l_1) + kg_vwakes_seg
                  end if
                  
                  if (node2 .ne. 0) then
                    indicesq_nodes_ij_l_2 = this%surfaces(this%wakes(k)%surface)%indicesq(this%surfaces(this%wakes(k)%surface)%nodes(node2)%indicesq)
            
                    kg_vwakes_seg = this%airflow%density * this%surfaces(i)%rings(j)%area * &
                        matmul(outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%deltav), kvwv_matrices%kvs3)                  
                  
                    k123_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),kg_vwakes_seg)
                    this%dfae(indicesq_rings_ij,indicesq_nodes_ij_l_2) = this%dfae(indicesq_rings_ij,indicesq_nodes_ij_l_2) + k123_temp
                    if (boolean_output_dfae_ring)  dfae_ring(indicesp_rings_ij,indicesq_nodes_ij_l_2) = dfae_ring(indicesp_rings_ij,indicesq_nodes_ij_l_2) + kg_vwakes_seg
                  end if
                end do
              end do
            end if
            !
            !  tangent matrix depending on geometry of ring j
            knn_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),kg_A + kg_n + kg_Dv_lvec + kg_vwakes)
            this%dfae(indicesq_rings_ij,indicesq_rings_ij) = this%dfae(indicesq_rings_ij,indicesq_rings_ij) + knn_temp

            !  tangent matrix depending on geometry of ring j
            if (boolean_output_dfae_ring)  dfae_ring(indicesp_rings_ij,indicesq_rings_ij) = dfae_ring(indicesp_rings_ij,indicesq_rings_ij) + kg_A + kg_n + kg_Dv_lvec + kg_vwakes

            !< linearization of vsurf
            kg_vsurf = - this%airflow%density * this%surfaces(i)%rings(j)%area * matmul(outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%deltav),this%surfaces(i)%rings(j)%kvsurf)

            !  tangent matrix depending on velocity of ring j
            knn_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),kg_vsurf)
            this%dfae(indicesq_rings_ij,indicesv_rings_ij) = this%dfae(indicesq_rings_ij,indicesv_rings_ij) + knn_temp

            !  tangent matrix depending on velocity of ring j
            if (boolean_output_dfae_ring) dfae_ring(indicesp_rings_ij,indicesv_rings_ij) = dfae_ring(indicesp_rings_ij,indicesv_rings_ij) + kg_vsurf

            !< calculate tangent matrix due to pressure differences
            do  k = 1, this%nsurfaces
              ! loop over all surface rings
              do l = 1, this%surfaces(k)%nrings
                indicesq_rings_kl = this%surfaces(k)%indicesq(this%surfaces(k)%rings(l)%indicesq)

                !< linearization of vsv
                call this%surfaces(k)%rings(l)%tangent_vortexringvelocity(kvs1, kvs2, this%surfaces(i)%rings(j)%center)

                !< linearization of vsv: derivative with respect to control point of ring j
                kg1 = this%airflow%density * this%surfaces(i)%rings(j)%area * this%surfaces(k)%rings(l)%circulation * &
                              matmul(outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%deltav), matmul(kvs1,this%surfaces(i)%rings(j)%Tcp))
                knn_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),kg1)
                this%dfae(indicesq_rings_ij,indicesq_rings_ij) = this%dfae(indicesq_rings_ij,indicesq_rings_ij) + knn_temp

                ! tangent matrix with respect to ring j
                if (boolean_output_dfae_ring) dfae_ring(indicesp_rings_ij,indicesq_rings_ij) = dfae_ring(indicesp_rings_ij,indicesq_rings_ij) + kg1

                !< linearization of vsv: derivative with respect to ring nodes of ring l
                kg1 = this%airflow%density * this%surfaces(i)%rings(j)%area * this%surfaces(k)%rings(l)%circulation * &
                              matmul(outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%deltav),kvs2)
                knn_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),kg1)
                this%dfae(indicesq_rings_ij,indicesq_rings_kl) = this%dfae(indicesq_rings_ij,indicesq_rings_kl) + knn_temp

                ! tangent matrix with respect to ring k
                if (boolean_output_dfae_ring) dfae_ring(indicesp_rings_ij,indicesq_rings_kl) = dfae_ring(indicesp_rings_ij,indicesq_rings_kl) + kg1
              end do

              ! loop over all surface nodes
              do l = 1, this%surfaces(k)%nnodes
                indicesq_nodes_kl = this%surfaces(k)%indicesq(this%surfaces(k)%nodes(l)%indicesq)
                indicesv_nodes_kl = this%tncoordinates + indicesq_nodes_kl

                !< linearization of vsv: derivative of ring circulation DG = dG/dq Dq
                kvsv_G_q = 0.0d0
                kvsv_G_v = 0.0d0
                do m = 1,this%tnrings
                    kvsv_G_q = kvsv_G_q + outer(this%surfaces(i)%rings(j)%avector(3*(m-1)+1:3*(m-1)+3),this%dG(m,indicesq_nodes_kl))
                    kvsv_G_v = kvsv_G_v + outer(this%surfaces(i)%rings(j)%avector(3*(m-1)+1:3*(m-1)+3),this%dG(m,indicesv_nodes_kl))
                end do
                k33_vsv_DG = this%airflow%density * this%surfaces(i)%rings(j)%area * outer(this%surfaces(i)%rings(j)%normal,matmul(this%surfaces(i)%rings(j)%deltav,kvsv_G_q))

                !< linearization of Dv: derivative of segment circulation DGamma = dGamma/dq Dq
                kdv_dcirc_q = 0.0d0
                kdv_dcirc_v = 0.0d0
                do m = 1, 4
                  kdv_dcirc_q = kdv_dcirc_q + outer( factor_tvec(3*(m-1)+1:3*(m-1)+3), this%surfaces(i)%DGamma(indicessegmentcirc(m), indicesq_nodes_kl))
                  kdv_dcirc_v = kdv_dcirc_v + outer( factor_tvec(3*(m-1)+1:3*(m-1)+3), this%surfaces(i)%DGamma(indicessegmentcirc(m), indicesv_nodes_kl))
                end do

                kdv_dcirc_q = - matmul(skew(this%surfaces(i)%rings(j)%normal), kdv_dcirc_q)/this%surfaces(i)%rings(j)%area
                k33_Dv_DG = this%airflow%density * this%surfaces(i)%rings(j)%area * &
                          outer(this%surfaces(i)%rings(j)%normal,matmul(this%surfaces(i)%rings(j)%vmean-this%surfaces(i)%rings(j)%vsurf,kdv_dcirc_q))

                !< linearization of G-G_old: DG = dG/dq Dq
                if (this%boolean_unsteady_term) then
                    k33_G_DG = - this%airflow%density * this%surfaces(i)%rings(j)%area / this%deltat * outer(this%surfaces(i)%rings(j)%normal,this%dG(i1+j,indicesq_nodes_kl))
                    if (this%tcounter .eq. 0) k33_G_DG = 0.0d0
                end if

                k123_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),k33_vsv_DG + k33_Dv_DG + k33_G_DG)
                this%dfae(indicesq_rings_ij,indicesq_nodes_kl) = this%dfae(indicesq_rings_ij,indicesq_nodes_kl) + k123_temp

                ! tangent matrix with respect to ring k
                if (boolean_output_dfae_ring) dfae_ring(indicesp_rings_ij,indicesq_nodes_kl) = dfae_ring(indicesp_rings_ij,indicesq_nodes_kl) + k33_vsv_DG + k33_Dv_DG + k33_G_DG

                !< linearization of vsv: derivative with respect to ring circulation DG = dG/dv dv
                k33_vsv_DG = this%airflow%density * this%surfaces(i)%rings(j)%area * &
                                matmul(outer(this%surfaces(i)%rings(j)%normal,this%surfaces(i)%rings(j)%deltav),kvsv_G_v)

                !< linearization of Dv: derivative with respect to segment circulation DGamma = dGamma/dv dv
                kdv_dcirc_v = - matmul(skew(this%surfaces(i)%rings(j)%normal), kdv_dcirc_v)/this%surfaces(i)%rings(j)%area
                k33_Dv_DG = this%airflow%density * this%surfaces(i)%rings(j)%area * &
                          outer(this%surfaces(i)%rings(j)%normal,matmul(this%surfaces(i)%rings(j)%vmean-this%surfaces(i)%rings(j)%vsurf,kdv_dcirc_v))

                !< linearization of G-G_old: DG = dG/dv Dv
                if (this%boolean_unsteady_term) then
                    k33_G_DG = - this%airflow%density * this%surfaces(i)%rings(j)%area / this%deltat * outer(this%surfaces(i)%rings(j)%normal,this%dG(i1+j,indicesv_nodes_kl))
                    if (this%tcounter .eq. 0) k33_G_DG = 0.0d0
                end if
                k123_temp = matmul(transpose(this%surfaces(i)%rings(j)%Tcp),k33_vsv_DG + k33_Dv_DG + k33_G_DG)
                this%dfae(indicesq_rings_ij,indicesv_nodes_kl) = this%dfae(indicesq_rings_ij,indicesv_nodes_kl) + k123_temp

                ! tangent matrix with respect to ring k
                if (boolean_output_dfae_ring) dfae_ring(indicesp_rings_ij,indicesv_nodes_kl) = dfae_ring(indicesp_rings_ij,indicesv_nodes_kl) + (k33_vsv_DG + k33_Dv_DG + k33_G_DG)

              end do
            end do

          end do

          i1 = i1 + this%surfaces(i)%nrings

      end do

      ! writing result matrices
      if (boolean_output_dfae_ring) then
        this%dfae_ring = dfae_ring
        call writeRealVectorToFile(this%fae,9999,'fae.dres')
        call writeRealMatrixToFile(this%dfae,9999,'dfae.dres')
        call writeRealMatrixToFile(this%dG,9999,'dG.dres')
        call writeRealMatrixToFile(dfae_ring,9999,'dfae_ring.dres')
      end if

    end if

    return
  end subroutine modeldfae

!< Subroutine for linearization of non-penetration condition
  subroutine linearized_non_penetration_condition(this)
    implicit none

    class(model_aero), intent(inout) :: this
    integer :: i, j, k, l, i1, info
    integer :: indicesq_ae_global_ij(12), indicesq_ae_global_kl(12), indicesv_ae_global_ij(12)
    integer :: size_wake_segments
    real(kind = 8) :: kvs1(3,3), kvs2(3,12), kvwv(3,3)
    real(kind = 8), allocatable :: amatrix(:,:)
    type(tangent_matrices3_3x3):: kvwv_matrices
    integer :: indicesq_nodes_ij_l_1(3), indicesq_nodes_ij_l_2(3)
    integer :: node1, node2, segment_ID

    allocate(amatrix(this%tnrings,this%tnrings))
    
    this%dAG  = 0.0d0
    this%drhs = 0.0d0
    this%dG   = 0.0d0
    amatrix   = 0.0d0

  !< calculate linearization of circulation on global level and of aerodynamical forces on ring level
    i1   = 0
    do i = 1, this%nsurfaces
        do j = 1, this%surfaces(i)%nrings
          indicesq_ae_global_ij = this%surfaces(i)%indicesq(this%surfaces(i)%rings(j)%indicesq)
          indicesv_ae_global_ij = this%tncoordinates + this%surfaces(i)%indicesq(this%surfaces(i)%rings(j)%indicesq)

      !< linearization of surface geometry: area, normal, segment tangent: can be done in model_aero%actualize
          call this%surfaces(i)%rings(j)%tangent_vortexringgeometry()

      !< linearization of vsurf
          this%surfaces(i)%rings(j)%kvsurf = this%surfaces(i)%rings(j)%Tcp

      !< linearization of vwakes with respect to bounded-vortex sheets i ring coordinates rj for all segments
          kvwv = 0.0d0
          !!$OMP parallel do reduction(+:kvwv) &
          !!$OMP private(k, l, size_wake_segments, kvwv_matrices)
          do k = 1, this%nwakes
            size_wake_segments = 0
            if (allocated(this%wakes(k)%segments_t)) size_wake_segments = size(this%wakes(k)%segments_t)
            do l = 1, size_wake_segments
              kvwv_matrices = this%wakes(k)%segments(this%wakes(k)%segments_t(l))%tangent_velocity(this%surfaces(i)%rings(j)%center)
              kvwv = kvwv + kvwv_matrices%kvs1
            end do
          end do
          !!$omp end parallel do
          this%surfaces(i)%rings(j)%kvwakes = matmul(kvwv,this%surfaces(i)%rings(j)%Tcp)

      !< linearization of right hand side: dbi = d (- (vinf + vwv - vsurf)*ni )
          !< area normal
          this%drhs(i1 + j, indicesq_ae_global_ij) = this%drhs(i1 + j, indicesq_ae_global_ij) + &
            - ( matmul( this%surfaces(i)%rings(j)%vinfinity + this%surfaces(i)%rings(j)%vwakes -  this%surfaces(i)%rings(j)%vsurf, this%surfaces(i)%rings(j)%knormal) )

          !< wakes with repsect to bounded-vortex sheets i ring j control point coordinate
          this%drhs(i1 + j, indicesq_ae_global_ij) = this%drhs(i1 + j, indicesq_ae_global_ij) + &
            - ( matmul( this%surfaces(i)%rings(j)%normal, this%surfaces(i)%rings(j)%kvwakes) )

          !< with respect to kinematic velocity
          this%drhs(i1 + j, indicesv_ae_global_ij) = this%drhs(i1 + j, indicesv_ae_global_ij) + &
                matmul( this%surfaces(i)%rings(j)%normal, this%surfaces(i)%rings(j)%kvsurf )

          !< linearization of segments coordinates; additional terms due to updated first wake row
          if (this%tcounter .ne. 0 .and. this%boolean_lin_wake_row .eqv. .TRUE.) then
            do k = 1, this%nwakes
              do l = 1, size(this%wakes(k)%segments_ID_first_row)
                segment_ID = this%wakes(k)%segments_ID_first_row(l)
                kvwv_matrices = this%wakes(k)%segments(segment_ID)%tangent_velocity(this%surfaces(i)%rings(j)%center)
              
                !< assining to global nodes of bounded-vortex sheet
                node1 = this%wakes(k)%segments(segment_ID)%connectivity_to_bounded_sheet(1)
                node2 = this%wakes(k)%segments(segment_ID)%connectivity_to_bounded_sheet(2)
                if (node1 .ne. 0) then
                  indicesq_nodes_ij_l_1 = this%surfaces(this%wakes(k)%surface)%indicesq(this%surfaces(this%wakes(k)%surface)%nodes(node1)%indicesq)
                  this%drhs(i1 + j, indicesq_nodes_ij_l_1) = this%drhs(i1 + j, indicesq_nodes_ij_l_1) + &
                      - ( matmul( this%surfaces(i)%rings(j)%normal, kvwv_matrices%kvs2) )
                end if
                if (node2 .ne. 0) then
                  indicesq_nodes_ij_l_2 = this%surfaces(this%wakes(k)%surface)%indicesq(this%surfaces(this%wakes(k)%surface)%nodes(node2)%indicesq)
                  this%drhs(i1 + j, indicesq_nodes_ij_l_2) = this%drhs(i1 + j, indicesq_nodes_ij_l_2) + &
                      - ( matmul( this%surfaces(i)%rings(j)%normal, kvwv_matrices%kvs3) )
                end if
              end do
            end do
          end if
          
      !< linearization of left hand-side in non-penetration condition d(A*G)
          this%dAG(i1 + j, indicesq_ae_global_ij) = this%dAG(i1 + j, indicesq_ae_global_ij) + &
                        matmul(this%surfaces(i)%rings(j)%vsurfaces,this%surfaces(i)%rings(j)%knormal)

          do  k = 1, this%nsurfaces
            do l = 1, this%surfaces(k)%nrings
              indicesq_ae_global_kl = this%surfaces(k)%indicesq(this%surfaces(k)%rings(l)%indicesq)

              call this%surfaces(k)%rings(l)%tangent_vortexringvelocity(kvs1, kvs2, this%surfaces(i)%rings(j)%center)

              this%dAG(i1 + j, indicesq_ae_global_ij) = this%dAG(i1 + j, indicesq_ae_global_ij) +  &
                    this%surfaces(k)%rings(l)%circulation * matmul(this%surfaces(i)%rings(j)%normal, matmul(kvs1,this%surfaces(i)%rings(j)%Tcp))

              this%dAG(i1 + j, indicesq_ae_global_kl) = this%dAG(i1 + j, indicesq_ae_global_kl) +  &
                    this%surfaces(k)%rings(l)%circulation * matmul(this%surfaces(i)%rings(j)%normal, kvs2)
            end do
          end do

        end do
        i1 = i1 + this%surfaces(i)%nrings
    end do

  !< solving for global dG = A^1*(drhs - dAG)
    this%dG = this%drhs - this%dAG
    amatrix = this%amatrix
    call dgesv(this%tnrings, 2*this%tncoordinates, amatrix, this%tnrings, this%ipiv, this%dG, this%tnrings, info)

  !< calculating Dgammai = Trsi*DGi for each surface
    do i = 1, this%nsurfaces
      this%surfaces(i)%DGamma = 0.0d0
      if (this%surfaces(i)%boolean_lifting_surface .eqv. .TRUE.) then
        !< extract only delta ring circulations for ring of surface i
        this%surfaces(i)%DG = this%dG(this%surfaces(i)%indicesc,1:2*this%tncoordinates)
        !< sparse multiplication to get delta segment circulation of surface i
        call mkl_dcsrmm('N', this%surfaces(i)%nsegments, 2*this%tncoordinates, this%surfaces(i)%nrings, 1.0d0, 'G  F  ', &
            this%surfaces(i)%sparse_Trs%csr_values, this%surfaces(i)%sparse_Trs%csr_columns, &
            this%surfaces(i)%sparse_Trs%csr_rowIndices(1:size(this%surfaces(i)%sparse_Trs%csr_rowIndices) - 1), &
            this%surfaces(i)%sparse_Trs%csr_rowIndices(2:size(this%surfaces(i)%sparse_Trs%csr_rowIndices)), &
            this%surfaces(i)%DG, this%surfaces(i)%nrings, 0.0d0, this%surfaces(i)%DGamma, this%surfaces(i)%nsegments)
      end if
    end do

    return
  end subroutine linearized_non_penetration_condition

!< Subroutine to calculate aerodynamic force vector and tangential matrix
  subroutine modelfae(this)
    implicit none

    class(model_aero), intent(inout) :: this
    integer :: i, j

    !< calculate aerodynamical force due to pressure
    this%fae = 0.0d0

    do i = 1, this%nsurfaces
      if (this%surfaces(i)%boolean_lifting_surface) then
        do j = 1, this%surfaces(i)%nrings
            call this%surfaces(i)%rings(j)%aerodynamicloads(this%circulations, this%surfaces(i)%segmentscirculations, this%surfaces(i)%oldsegmentscirculations, this%tnrings, this%deltat, this%tcounter, this%airflow%density, this%boolean_unsteady_term)
            this%fae(this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp)) = this%surfaces(i)%rings(j)%fae
        end do
      end if
    end do

    ! transformation from control points to aero nodes
    call mkl_dcsrmv('T', 3*this%tnrings, this%tncoordinates, 1.0d0, 'G  F  ', &
    this%sparse_Tcp3%csr_values, this%sparse_Tcp3%csr_columns, &
    this%sparse_Tcp3%csr_rowIndices(1:size(this%sparse_Tcp3%csr_rowIndices) - 1), &
    this%sparse_Tcp3%csr_rowIndices(2:size(this%sparse_Tcp3%csr_rowIndices)), &
    this%fae, 0.0d0, this%faen)
    
    return
  end subroutine modelfae

!< Function to determine amplitude of a function
  function evaluateamplitude(this,time) result(amplitude)
    implicit none

    class(model_aero), intent(inout) :: this

    real(kind = 8) :: time
    real(kind = 8) :: amplitude

    select case (this%airflow%sort)
      !< dirac delta. Only works for exact times
      case('delta')
         if ((time >= this%airflow%duration) .AND. (time <= this%airflow%duration + 1.0d0)) then
            amplitude = 1.0d0
         else
            amplitude = 0.0d0
         end if

      ! step function, normalized
      case('heaviside')
        if( time > this%airflow%duration) then
          amplitude = 1.0d0
        else
          amplitude = 0.0d0
        endif

      !< step function, normalized
      case('constant')
        if( time >= 0.0d0) then
          amplitude = 1.0d0
        else
          amplitude = 0.0d0
        endif

      !< linear increase
      case('linear')
        amplitude = time/this%airflow%duration

      !< linear increase until duration, then amplitude = 0.0d0
      case('linear0')
        if (time <= this%airflow%duration) then
          amplitude = time/this%airflow%duration
        else
          amplitude = 0.0d0
        end if

      !< linear increase until duration and then constant
      case('linearC')
        if (time <= this%airflow%duration) then
          amplitude = time/this%airflow%duration
        else
          amplitude = 1.0d0
        end if

      !< linear increase until half of the duration, then linear decay
      case('triangle')
        if (time < this%airflow%duration/2.0d0) then
          amplitude = time/(this%airflow%duration/2.0d0)
        else if ((time>=this%airflow%duration/2.0d0) .and. (time<this%airflow%duration)) then
          amplitude = (1.0d0-(time-this%airflow%duration/2.0d0)/(this%airflow%duration/2.0d0))
        else
          amplitude = 0.0d0
        end if

      !< constant -> linear increase -> 1.3*constant
      case('event')
        if (time < (0.5d0*this%airflow%duration)) then
          amplitude = 1.0d0
        else if (time >= (0.5d0*this%airflow%duration) .and. (time < (0.55d0*this%airflow%duration))) then
          amplitude = (6.0d0/(this%airflow%duration)) * time - 2.0d0
        else
          amplitude = 1.3d0
        end if  
        
      !< constant -> linear decrease -> 0.7*constant
      case('event_low')
        if (time < (0.5d0*this%airflow%duration)) then
          amplitude = 1.0d0
        else if (time >= (0.5d0*this%airflow%duration) .and. (time < (0.55d0*this%airflow%duration))) then
          amplitude = (-6.0d0/(this%airflow%duration)) * time + 4.0d0
        else
          amplitude = 0.7d0
        end if
        
        
      !< inverted heaviside
      case('square')
        if (time <= this%airflow%duration) then
          amplitude = 1.0d0
        else if (time > this%airflow%duration) then
          amplitude = 0.0d0
        end if
      !< default case
      case default
        amplitude = 0.0d0
        print*, 'Wanring: wind field function/type not specified!'
        this%boolean_abort_airflow = .TRUE.
    end select

    !< scale normalized load amplitude by intensity
    amplitude = amplitude *this%airflow%intensity

    return
  end function evaluateamplitude


!< Subroutine to read kinematic of bounded vortex sheeet from file
  subroutine read_kinematic(this)
    implicit none
    class(model_aero), intent(inout) :: this
    integer :: i, io_errorFile, ireason, UnIn
    integer :: nsteps
    character(:), allocatable :: input_filename

    !< reading solution file from DeSiO-Structure solution q
    input_filename = trim(this%output_filename)
    call GetNewUnit (UnIn)
    open(unit = UnIn, file = input_filename // '_uvlm_qs.dres', status = 'old', iostat = io_errorFile)
    !< Determining number of line in input file
    if (io_errorFile .eq. 0) then
      print*,'loading kinematic coordinate file...', input_filename // '_uvlm_qs.dres'
      this%bool_const_a_matrix = .FALSE.
      nsteps = 0
      read(UnIn, *, iostat = ireason)
      do i = 1, this%nsteps
          if (ireason .lt. 0) exit
          read(UnIn, *, iostat = ireason)
          nsteps = nsteps + 1
      end do
      close(unit = UnIn)

      allocate(this%arrqs_t(this%nsteps,this%tncoordinates))
      this%arrqs_t(:,:) = 0.0d0
      call GetNewUnit (UnIn)
      open(unit = UnIn, file = input_filename // '_uvlm_qs.dres', status = 'old', iostat = io_errorFile)
      do i = 1, this%nsteps
        if (i .gt. nsteps) then
          this%arrqs_t(i,:) = this%arrqs_t(nsteps,:)
        else
          read(UnIn, *) this%arrqs_t(i,:)
        end if
      end do
      close(unit = UnIn)
    else
      print*,'no kinematic coordinate file available...'
      print*,'no kinematic velocity file available...'
      print*,'using constant coefficient matrix A...'
      close(unit = UnIn)
      return
    end if

    !< read velocity at surface control point from DeSiO-Structure solution v
    call GetNewUnit (UnIn)
    open(unit = UnIn, file = input_filename // '_uvlm_vs.dres', status = 'old', iostat = io_errorFile)
    if (io_errorFile .eq. 0) then
      print*,'loading kinematic velocity file...', input_filename // '_uvlm_vs.dres'
      allocate(this%arrvs_t(this%nsteps,this%tnvelocities))
      this%arrvs_t(:,:) = 0.0d0
      do i = 1, this%nsteps
        if (i .gt. nsteps) then
          this%arrvs_t(i,:) = this%arrvs_t(nsteps,:)
        else
          read(UnIn, *) this%arrvs_t(i,:)
        end if
      end do
    end if
    close(unit = UnIn)

    return
  end subroutine read_kinematic

!< Subroutine for opening result files
  subroutine output_open(this)
    implicit none

    class(model_aero), intent(in) :: this

    ! integer(8), parameter :: block_size = 2 ** 200
    integer :: j, i

    character(:), allocatable :: output_filename, char_i
    character(len=1024) :: char_t
    character(len=7) :: str_status = 'replace'
    character(len=6) :: str_position = 'rewind'

    1000 format(1000000f30.5)
    allocate(output_filename , source =  trim(adjustl(this%output_filename)))

    !< writing surface informations
    open(unit = 5001, file = output_filename // '_uvlm_models.dres', action = 'write', status = 'replace')!,  BUFFERED='YES', BLOCKSIZE=block_size)
        write(5001, *) this%nsurfaces
        do i = 1, this%nsurfaces
          write(5001, *) this%surfaces(i)%nnodes, this%nsteps, this%nsteps, this%surfaces(i)%nrings, this%surfaces(i)%nnodes
          do j = 1, this%surfaces(i)%nrings
            write(5001, *) this%surfaces(i)%rings(j)%connectivity
          end do
        end do
    close(5001)

    if (this%boolean_oldsim) then
      if (this%output_filename .eq. this%oldsimfilename) then
        str_status   = 'unknown'
        str_position = 'append'
      end if
    end if

    open(unit = 5002, file = output_filename // '_uvlm_circs_cp.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    open(unit = 5003, file = output_filename // '_uvlm_vs_nodal.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    open(unit = 5004, file = output_filename // '_uvlm_dps_cp.dres',   action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    open(unit = 5005, file = output_filename // '_uvlm_qs_nodal.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    open(unit = 5006, file = output_filename // '_uvlm_circs_nodal.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    open(unit = 5007, file = output_filename // '_uvlm_dps_nodal.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    open(unit = 5008, file = output_filename // '_uvlm_areanormal_nodal.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    !open(unit = 5009, file = output_filename // '_uvlm_vfluid_nodal.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)

    !open(unit = 5101, file = output_filename // '_uvlm_dv_cp.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    open(unit = 5102, file = output_filename // '_uvlm_vm_cp.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    !open(unit = 5103, file = output_filename // '_uvlm_vs_cp.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    !open(unit = 5104, file = output_filename // '_uvlm_vsv_cp.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)
    !open(unit = 5105, file = output_filename // '_uvlm_vwv_cp.dres', action = 'write', status = str_status, POSITION = str_position)!, BUFFERED='YES', BLOCKSIZE=block_size)

    !< checking, if old wake information file is available or not
    open(unit = 6001, file = output_filename // '_uvlm_modelw.dres', action = 'write', status = 'replace')!, BUFFERED='YES', BLOCKSIZE=block_size)
      write(6001, *) this%nwakes
      do i = 1, this%nwakes
        write(6001, *) this%wakes(i)%nnodess, this%wakes(i)%nstepstotal, this%wakes(i)%nrowstoconsider, this%wakes(i)%nrings, this%wakes(i)%nnodes
        do j = 1, this%wakes(i)%nrings
          write(6001, *) this%wakes(i)%rings(j)%connectivity
        end do
      end do
    close(6001)

    do i = 1,this%nwakes
      write(char_t, '(i5)')  i
      allocate(char_i , source =  trim(adjustl(char_t)))
      open(unit = 6300+i, file = output_filename // '_uvlm_qw' // char_i // '.dres', action = 'write', status = str_status, POSITION = str_position)
      open(unit = 6400+i, file = output_filename // '_uvlm_vw' // char_i // '.dres', action = 'write', status = str_status, POSITION = str_position)
      open(unit = 6500+i, file = output_filename // '_uvlm_circw' // char_i // '.dres', action = 'write', status = str_status, POSITION = str_position)
      deallocate(char_i)
    end do

    !< writing time information
    open(unit = 7001, file = output_filename // '_uvlm_t.dres', action = 'write', status = str_status, POSITION = str_position)

    return
  end subroutine output_open

!< Subroutine for writing result files
  subroutine output_write(this)
    implicit none

    class(model_aero), intent(inout) :: this
    integer :: i, j

1000 format(1000000f30.15)
1100 format(1000000E29.20E3)

    !< global results at control points
    do i = 1, this%nsurfaces
      do j = 1,this%surfaces(i)%nrings
        this%deltaps(this%surfaces(i)%indicesc(j))                                     = this%surfaces(i)%rings(j)%deltaps
        this%areanormal(this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp)) = this%surfaces(i)%rings(j)%normal

        !< velocities at control points
        this%deltav(this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp)) = this%surfaces(i)%rings(j)%deltav
        this%vmean(this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp))  = this%surfaces(i)%rings(j)%vmean
        this%vs(this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp))     = this%surfaces(i)%rings(j)%vsurf
        this%vsv(this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp))    = this%surfaces(i)%rings(j)%vsurfaces
        this%vwv(this%surfaces(i)%indicesp(this%surfaces(i)%rings(j)%indicesp))    = this%surfaces(i)%rings(j)%vwakes

      end do
    end do

    !< transfering results from control point to nodes
    call mkl_dcsrmv('T', this%tnrings, this%tnnodes, 1.0d0, 'G  F  ', &
    this%sparse_Tcp1%csr_values, this%sparse_Tcp1%csr_columns, &
    this%sparse_Tcp1%csr_rowIndices(1:size(this%sparse_Tcp1%csr_rowIndices) - 1), &
    this%sparse_Tcp1%csr_rowIndices(2:size(this%sparse_Tcp1%csr_rowIndices)), &
    this%deltaps, 0.0d0, this%deltaps_nodal)

    call mkl_dcsrmv('T', this%tnrings, this%tnnodes, 1.0d0, 'G  F  ', &
    this%sparse_Tcp1%csr_values, this%sparse_Tcp1%csr_columns, &
    this%sparse_Tcp1%csr_rowIndices(1:size(this%sparse_Tcp1%csr_rowIndices) - 1), &
    this%sparse_Tcp1%csr_rowIndices(2:size(this%sparse_Tcp1%csr_rowIndices)), &
    this%circulations, 0.0d0, this%circulations_nodal)

    call mkl_dcsrmv('T', 3*this%tnrings, this%tncoordinates, 1.0d0, 'G  F  ', &
    this%sparse_Tcp3%csr_values, this%sparse_Tcp3%csr_columns, &
    this%sparse_Tcp3%csr_rowIndices(1:size(this%sparse_Tcp3%csr_rowIndices) - 1), &
    this%sparse_Tcp3%csr_rowIndices(2:size(this%sparse_Tcp3%csr_rowIndices)), &
    this%areanormal, 0.0d0, this%areanormal_nodal)

    call mkl_dcsrmv('T', 3*this%tnrings, this%tnvelocities, 1.0d0, 'G  F  ', &
    this%sparse_Tcp3%csr_values, this%sparse_Tcp3%csr_columns, &
    this%sparse_Tcp3%csr_rowIndices(1:size(this%sparse_Tcp3%csr_rowIndices) - 1), &
    this%sparse_Tcp3%csr_rowIndices(2:size(this%sparse_Tcp3%csr_rowIndices)), &
    this%vmean, 0.0d0, this%vmean_nodal)

    !< writing bounded vortex sheet (vortex ring, vortex lattice) informations
    write(5002, 1000) this%circulations(:)
    write(5004, 1000) this%deltaps(:)

    write(5003, 1000) this%vs_t(:)
    write(5005, 1000) this%qs_t(:)
    write(5006, 1000) this%circulations_nodal(:)
    write(5007, 1000) this%deltaps_nodal(:)
    write(5008, 1000) this%areanormal_nodal(:)
    !write(5009, 1000) this%vmean_nodal(:)

    !< velocities at control point
    !write(5101, 1000) this%deltav(:)
    write(5102, 1000) this%vmean(:)
    !write(5103, 1000) this%vs(:)
    !write(5104, 1000) this%vsv(:)
    !write(5105, 1000) this%vwv(:)

    !< writing wake informations (unbounded vortex lattice)
    do i = 1,this%nwakes
      write(6300+i, 1000) this%wakes(i)%q_t(1:3*this%wakes(i)%nnodess*(this%wakes(i)%tcounter+1))

      rewind(6400+i) ! jump to initial position of the file
      write(6400+i, 1000) this%wakes(i)%v_t(1:3*this%wakes(i)%nnodess*(this%wakes(i)%tcounter+1))

      rewind(6500+i) ! jump to initial position of the file
      write(6500+i, 1000) this%wakes(i)%ringscirculations(1:this%wakes(i)%nringss*(this%wakes(i)%tcounter))

    end do

    !< write time informations
    write(7001, 1000) this%time, real(this%tcounter), this%delta_sim_time

    return
  end subroutine output_write

!< Subroutine for writing result files
  subroutine output_close(this)
    implicit none

    class(model_aero), intent(in) :: this
    integer :: i

    close(5002)
    close(5003)
    close(5004)
    close(5005)
    close(5006)
    close(5007)
    close(5008)
    !close(5009)

    !close(5101)
    close(5102)
    !close(5103)
    !close(5104)
    !close(5105)

    do i = 1,this%nwakes
      close(6300+i)
      close(6400+i)
      close(6500+i)
    end do

    close(7001)
    return
  end subroutine output_close

end module class_model_aero
