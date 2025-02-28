module class_model_fsi
  use class_model_aero, only: model_aero
  use class_model_structure, only: model_structure
  use my_constants_aero, only : pi, i_33
  use my_constants_structure, only: e1, e2, e3
  use my_math_aero, only: cross, skew, outer, eye
  use my_fsi_functions
  use sparse_function, only: single_sparse_matrix, create_rowcol_format, sparse_initialize
  implicit none
  
  type :: tfsi
      integer, allocatable :: structural_body(:), structural_nodes(:)
      integer :: fluid_surface
      real(kind = 8) :: search_radius
      real(kind = 8), allocatable :: arr_searchradius(:,:)
      !character(5) :: strtype, lfsort
      character(10) :: strtype, lfsort
      character(:), allocatable ::str_searchradius_filename
  end type tfsi
  
  type :: bounded_surfaces
      real(kind = 8), allocatable :: DGamma(:,:), DGamma_dense(:,:), DG(:,:)
      integer, allocatable :: indices_qv_aero(:)
  end type
  
  type :: model_fsi
    integer :: nfsi
    type(tfsi), allocatable :: fsi(:)
    type(single_sparse_matrix) :: sparse_Tn, sparse_T2n
    type(bounded_surfaces), allocatable :: surfaces(:)
    logical :: boolean_abort = .FALSE.
    logical :: boolean_strongfsi = .TRUE.
    real(kind = 8), allocatable :: temp_fn(:), fae(:), fae_x(:), fae_w(:)
    real(kind = 8), allocatable :: dAG(:,:), drhs(:,:), dG(:,:)
    real(kind = 8), allocatable :: an(:)
    real(kind = 8) :: timef
    integer :: flag_dfae
    character(:), allocatable :: output_filename, lfsort, str_fsi_type

    !< for material aerodynamic loads 
    real( kind = 8), allocatable :: ml_aero_t(:)
    
    contains
      procedure :: constructor
      procedure :: actualize
      procedure :: readinginputs
      procedure :: calc_load_factor
      procedure :: amplitude_material_loads
    end type model_fsi
 
  contains

!< Subroutine for construct the model  
  subroutine constructor(this,this_structure,this_aero)
    use sparse_function
    implicit none

    class(model_fsi), intent(inout) :: this
    class(model_structure), intent(in) :: this_structure
    class(model_aero), intent(inout) :: this_aero
    integer :: i, j, k, l, nj, k1
    integer :: indicesq_st_global(12), indicesv_st_global(12), indicesq_ae_global(3), indicesp_ae_global(3), indices_temp(12)
    integer, allocatable :: indicesq_st_local(:)
    integer, allocatable :: arr_nodes(:), inz(:)
    real(kind = 8), dimension(3) :: xae, x0, d10, d20, d30
    real(kind = 8), allocatable :: arr_len(:), arr_xhi10(:), arr_xhi20(:), arr_xhi30(:), arr_weights(:), arr_d_ae_st(:,:), arr_d0(:,:), arr_a_ae_st(:,:)
    real(kind = 8) :: L_ae_st, search_radius
    real(kind = 8) :: d0(3), d_ae_st(3), n_ae_st(3), R_ae_st(3,3), a_ae_st(3)
            
    integer :: nn_st
    integer, allocatable :: arr_nodes_st(:)
    integer :: convert_job_coo_csr(8), info, nNonZeros
    integer :: r1
    integer :: nnodes, nnodes_st
    integer :: nodeIDgl, nfsi_strb

    convert_job_coo_csr    = 0 ! the matrix in the CSR format is converted to the coordinate format;
    convert_job_coo_csr(1) = 2 ! the matrix in the coordinate format is converted to the CSR format, and the column indices in CSR representation are sorted in the increasing order within each row.
    convert_job_coo_csr(2) = 1 ! one-based indexing for the matrix in CSR format is used.
    convert_job_coo_csr(3) = 1 ! one-based indexing for the matrix in coordinate format is used.
    convert_job_coo_csr(6) = 0 ! all arrays acsr, ja, ia are filled in for the output storage.
    
    this%nfsi = 0
    if (allocated(this%fsi)) then
      this%nfsi = size(this%fsi)
    end if
    
    if (this%nfsi .eq. 0) this%boolean_abort = .TRUE.
    
    allocate(this%temp_fn(this_aero%tncoordinates))
    this%temp_fn = 0.0d0

    allocate(this%fae(this_structure%tncoordinates))
    this%fae = 0.0d0
    
    allocate(this%fae_x(this_structure%tncoordinates))
    this%fae_x = 0.0d0

    allocate(this%fae_w(this_structure%tncoordinates))
    this%fae_w = 0.0d0

    allocate(this%an(this_aero%tncoordinates))
    this%an  = 0.0d0

    allocate(this%ml_aero_t(6*this_structure%nnodes12+3*this_structure%nnodes6))
    this%ml_aero_t(:) = 0.0d0

    ! initialize sparse matrices
    call sparse_initialize(this%sparse_Tn, this_aero%tncoordinates, this_structure%tncoordinates)
    call sparse_initialize(this%sparse_T2n, 2*this_aero%tncoordinates, 2*this_structure%tncoordinates)
    
    !< loop over all FSI (beam - surface)
    print*, 'Calculating sparse transformation matrix for FSI'
    do k = 1, this%nfsi
      
      if (this%fsi(k)%strtype .eq. 'beam') then
        
        !< loop over all nodes of surface
        do i = 1, this_aero%surfaces(this%fsi(k)%fluid_surface)%nnodes
          
          !< serach radius for the node
          search_radius = this%fsi(k)%arr_searchradius(i,2)
          
          !< getting aero node coordinates and global indices
          xae = this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%q_t
          indicesq_ae_global = this_aero%surfaces(this%fsi(k)%fluid_surface)%indicesq(this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%indicesq)

          nfsi_strb = 0
          if (allocated(this%fsi(k)%structural_body)) then
            nfsi_strb = size(this%fsi(k)%structural_body)
          end if

          !< loop over all nodes of beam structure
          do k1 = 1, nfsi_strb
            do j = 1, this_structure%beams(this%fsi(k)%structural_body(k1))%nnodes
              !< global node ID
              nodeIDgl = this_structure%beams(this%fsi(k)%structural_body(k1))%nodes(j)%node_structure%globalid
              
              !< getting structural node coordinates
              x0  = this_structure%q_0(this_structure%beams(this%fsi(k)%structural_body(k1))%indicesq12(12*(j-1)+1 :12*(j-1)+3))
              d10 = this_structure%q_0(this_structure%beams(this%fsi(k)%structural_body(k1))%indicesq12(12*(j-1)+4 :12*(j-1)+6))
              d20 = this_structure%q_0(this_structure%beams(this%fsi(k)%structural_body(k1))%indicesq12(12*(j-1)+7 :12*(j-1)+9))
              d30 = this_structure%q_0(this_structure%beams(this%fsi(k)%structural_body(k1))%indicesq12(12*(j-1)+10:12*(j-1)+12))
          
              !< length between node ae and node st
              L_ae_st = norm2(xae - x0)
          
              !< store only nodes inside search range
              if (L_ae_st .lt. search_radius) then
                call append_int_vec(arr_nodes,nodeIDgl,0)
                call append_real_vec(arr_len,L_ae_st)
                call append_real_vec(arr_xhi10,dot_product((xae - x0),d10))
                call append_real_vec(arr_xhi20,dot_product((xae - x0),d20))
                call append_real_vec(arr_xhi30,dot_product((xae - x0),d30))
                call append_real_vec(arr_weights,weight_factor(L_ae_st,search_radius,1))
              end if
            end do
          end do
          
          !< loop over all found beam nodes and fill sparse transformation matrix for global structure nodes to aero surface nodes
          nnodes = 0
          if (allocated(arr_nodes)) then
            nnodes = size(arr_nodes)
          else
            print*, 'Warning in DeSiO-FSI: no nodes coupled for aero node'
          end if
          
          do j = 1,nnodes
            indicesq_st_global = this_structure%nodes12(arr_nodes(j))%coordinates
            
            !< delta_x * ( f(q) * DX + f(q) * DV)
            !< for Tn
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
            ! for T2n = [Tn Tn]
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
          end do
        
          !< deallocating
          deallocate(arr_xhi10)
          deallocate(arr_xhi20)
          deallocate(arr_xhi30)
          deallocate(arr_nodes)
          deallocate(arr_weights)
          deallocate(arr_len)
        end do
      
      elseif (this%fsi(k)%strtype .eq. 'rigid_body') then
        
        !< loop over all nodes of surface
        do i = 1, this_aero%surfaces(this%fsi(k)%fluid_surface)%nnodes
          
          !< serach radius for the node
          search_radius = this%fsi(k)%arr_searchradius(i,2)
          
          !< getting aero node coordinates and global indices
          xae = this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%q_t
          indicesq_ae_global = this_aero%surfaces(this%fsi(k)%fluid_surface)%indicesq(this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%indicesq)

          nfsi_strb = 0
          if (allocated(this%fsi(k)%structural_body)) then
            nfsi_strb = size(this%fsi(k)%structural_body)
          end if

          !< loop over all nodes of rigid body structure
          do k1 = 1, nfsi_strb
            !< global node ID
            nodeIDgl = this_structure%bodies(this%fsi(k)%structural_body(k1))%node%node_structure%globalid
              
            !< getting structural node coordinates
            x0  = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(1 :3))
            d10 = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(4 :6))
            d20 = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(7 :9))
            d30 = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(10:12))
          
            !< length between node ae and node st
            L_ae_st = norm2(xae - x0)
          
            !< store only nodes inside search range
            if (L_ae_st .lt. search_radius) then
              call append_int_vec(arr_nodes,nodeIDgl,0)
              call append_real_vec(arr_len,L_ae_st)
              call append_real_vec(arr_xhi10,dot_product((xae - x0),d10))
              call append_real_vec(arr_xhi20,dot_product((xae - x0),d20))
              call append_real_vec(arr_xhi30,dot_product((xae - x0),d30))
              call append_real_vec(arr_weights,weight_factor(L_ae_st,search_radius,1))
            end if
          end do
          
          !< loop over all found beam nodes and fill sparse transformation matrix for global structure nodes to aero surface nodes
          nnodes = 0
          if (allocated(arr_nodes)) then
            nnodes = size(arr_nodes)
          else
            print*, 'Warning in DeSiO-FSI: no nodes coupled for aero node'
          end if
          
          do j = 1,nnodes
            indicesq_st_global = this_structure%nodes12(arr_nodes(j))%coordinates
            
            !< delta_x * ( f(q) * DX + f(q) * DV)
            !< for Tn
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
            ! for T2n = [Tn Tn]
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
          end do
        
          !< deallocating
          deallocate(arr_xhi10)
          deallocate(arr_xhi20)
          deallocate(arr_xhi30)
          deallocate(arr_nodes)
          deallocate(arr_weights)
          deallocate(arr_len)
        end do
      
      elseif (this%fsi(k)%strtype .eq. 'rigid_body') then
        
        !< loop over all nodes of surface
        do i = 1, this_aero%surfaces(this%fsi(k)%fluid_surface)%nnodes
          
          !< serach radius for the node
          search_radius = this%fsi(k)%arr_searchradius(i,2)
          
          !< getting aero node coordinates and global indices
          xae = this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%q_t
          indicesq_ae_global = this_aero%surfaces(this%fsi(k)%fluid_surface)%indicesq(this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%indicesq)

          nfsi_strb = 0
          if (allocated(this%fsi(k)%structural_body)) then
            nfsi_strb = size(this%fsi(k)%structural_body)
          end if

          !< loop over all nodes of beam structure
          do k1 = 1, nfsi_strb
            !< global node ID
            nodeIDgl = this_structure%bodies(this%fsi(k)%structural_body(k1))%node%node_structure%globalid
              
            !< getting structural node coordinates
            x0  = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(1:3))
            d10 = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(4:6))
            d20 = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(7:9))
            d30 = this_structure%q_0(this_structure%bodies(this%fsi(k)%structural_body(k1))%indicesq12(10:12))
          
            !< length between node ae and node st
            L_ae_st = norm2(xae - x0)
          
            !< store only nodes inside search range
            if (L_ae_st .lt. search_radius) then
              call append_int_vec(arr_nodes,nodeIDgl,0)
              call append_real_vec(arr_len,L_ae_st)
              call append_real_vec(arr_xhi10,dot_product((xae - x0),d10))
              call append_real_vec(arr_xhi20,dot_product((xae - x0),d20))
              call append_real_vec(arr_xhi30,dot_product((xae - x0),d30))
              call append_real_vec(arr_weights,weight_factor(L_ae_st,search_radius,1))
            end if
          end do
          
          !< loop over all found beam nodes and fill sparse transformation matrix for global structure nodes to aero surface nodes
          nnodes = 0
          if (allocated(arr_nodes)) then
            nnodes = size(arr_nodes)
          else
            print*, 'Warning in DeSiO-FSI: no nodes coupled for aero node'
          end if
          
          do j = 1,nnodes
            indicesq_st_global = this_structure%nodes12(arr_nodes(j))%coordinates
            
            !< for Tn
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
            ! for T2n
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(4:6),   i_33*arr_xhi10(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(7:9),   i_33*arr_xhi20(j)*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(10:12), i_33*arr_xhi30(j)*arr_weights(j)/sum(arr_weights))
            
          end do
        
          !< deallocating
          deallocate(arr_xhi10)
          deallocate(arr_xhi20)
          deallocate(arr_xhi30)
          deallocate(arr_nodes)
          deallocate(arr_weights)
          deallocate(arr_len)
        end do
        
      elseif (this%fsi(k)%strtype .eq. 'shell') then
      
        !< loop over all nodes of surface
        do i = 1, this_aero%surfaces(this%fsi(k)%fluid_surface)%nnodes

          !< serach radius for the node
          search_radius = this%fsi(k)%arr_searchradius(i,2)
          
          !< getting aero node coordinates and global indices
          xae = this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%q_t
          indicesq_ae_global = this_aero%surfaces(this%fsi(k)%fluid_surface)%indicesq(this_aero%surfaces(this%fsi(k)%fluid_surface)%nodes(i)%indicesq)
        
          nfsi_strb = 0
          if (allocated(this%fsi(k)%structural_body)) then
            nfsi_strb = size(this%fsi(k)%structural_body)
          end if

          !< loop over all nodes of shell structure
          do k1 = 1,nfsi_strb
            do j = 1, this_structure%shells(this%fsi(k)%structural_body(k1))%nnodes
              !< global node ID
              nodeIDgl = this_structure%shells(this%fsi(k)%structural_body(k1))%nodes(j)%node_structure%globalid

              !< getting structural node coordinates
              x0 = this_structure%q_0(this_structure%shells(this%fsi(k)%structural_body(k1))%indicesq6(6*(j-1)+1:6*(j-1)+3))
              d0 = this_structure%q_0(this_structure%shells(this%fsi(k)%structural_body(k1))%indicesq6(6*(j-1)+4:6*(j-1)+6))
          
              !< length and direction between node ae and node st
              a_ae_st = xae - x0
              L_ae_st = norm2(xae - x0)
                        
              !< store only nodes inside search range
              if (L_ae_st .lt. search_radius) then
                !d_ae_st = d0
                !if (L_ae_st .ne. 0.0d0) d_ae_st = (xae - x0)/L_ae_st 
                call append_int_vec(arr_nodes,j,0)
                call append_real_vec(arr_weights,weight_factor(L_ae_st,search_radius,1))
                call append_real_vec3(arr_a_ae_st,a_ae_st)
                !call append_real_vec(arr_len,L_ae_st)
                !call append_real_vec3(arr_d0,d0)
                !call append_real_vec3(arr_d_ae_st,d_ae_st)
              end if
            end do
          end do
          
          !< loop over all found beam nodes and fill sparse transformation matrix for global structure nodes to aero surface nodes
          nnodes = 0
          if (allocated(arr_nodes)) then
            nnodes = size(arr_nodes)
          else
            print*, 'Warning in DeSiO-FSI: no nodes coupled for aero node'
          end if
        
          do j = 1,nnodes
            !< calculating rotation matrix
            !d_ae_st = arr_d_ae_st(j,1:3)
            !d0      = arr_d0(j,1:3)
            !n_ae_st = cross(d0,d_ae_st)
            !R_ae_st = eye(3) + norm2(n_ae_st)*skew(n_ae_st) + (1.0d0 - dot_product(d0,d_ae_st))*matmul(skew(n_ae_st),skew(n_ae_st)) 
            indicesq_st_global(1:6)     = this_structure%nodes6(arr_nodes(j))%coordinates
            a_ae_st                     = arr_a_ae_st(j,1:3)
            this%an(indicesq_ae_global) = this%an(indicesq_ae_global) + arr_weights(j)/sum(arr_weights)*a_ae_st
            call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, indicesq_ae_global, indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            call create_rowcol_format(this%sparse_T2n, this_aero%tncoordinates + indicesq_ae_global, this_structure%tncoordinates + indicesq_st_global(1:3),   i_33*arr_weights(j)/sum(arr_weights))
            !call create_rowcol_format(this%sparse_Tn, indicesq_ae_global, indicesq_st_global(4:6), R_ae_st*arr_len(j)*arr_weights(j)/sum(arr_weights))
            !this%Tn(indicesq_ae_global,indicesq_st_global(1:3)) = i_33*arr_weights(j)/sum(arr_weights)
            !this%Tn(indicesq_ae_global,indicesq_st_global(4:6)) = R_ae_st*arr_len(j)*arr_weights(j)/sum(arr_weights)
          end do

          !< deallocating
          deallocate(arr_nodes)
          deallocate(arr_weights)
          deallocate(arr_a_ae_st)
          !deallocate(arr_len)
          !deallocate(arr_d_ae_st)
          !deallocate(arr_d0)
        end do
      
      end if
      
    end do
    
    !< convert transformation matrix Tn to csr sparse format
    allocate(this%sparse_Tn%csr_values(this%sparse_Tn%ncoords))
    allocate(this%sparse_Tn%csr_columns(this%sparse_Tn%ncoords))
    allocate(this%sparse_Tn%csr_rowIndices(this_aero%tncoordinates+1))
    this%sparse_Tn%csr_values     = 0.0d0
    this%sparse_Tn%csr_columns    = 0
    this%sparse_Tn%csr_rowIndices = 0
    !< convert coo format to csr format: call mkl_scsrcoo(job, n, acsr, ja, ia, nnz, acoo, rowind, colind, info)
    call mkl_dcsrcoo(convert_job_coo_csr, this_aero%tncoordinates, this%sparse_Tn%csr_values, &
                  this%sparse_Tn%csr_columns, this%sparse_Tn%csr_rowIndices, this%sparse_Tn%ncoords, & 
                  this%sparse_Tn%coo_values(1:this%sparse_Tn%ncoords), this%sparse_Tn%coo_rows(1:this%sparse_Tn%ncoords), this%sparse_Tn%coo_cols(1:this%sparse_Tn%ncoords+1), info)
    
    !< convert transformation matrix T2n to csr sparse format
    allocate(this%sparse_T2n%csr_values(this%sparse_T2n%ncoords))
    allocate(this%sparse_T2n%csr_columns(this%sparse_T2n%ncoords))
    allocate(this%sparse_T2n%csr_rowIndices(2*this_aero%tncoordinates+1))
    this%sparse_T2n%csr_values     = 0.0d0
    this%sparse_T2n%csr_columns    = 0
    this%sparse_T2n%csr_rowIndices = 0
    !< convert coo format to csr format: call mkl_scsrcoo(job, n, acsr, ja, ia, nnz, acoo, rowind, colind, info)
    call mkl_dcsrcoo(convert_job_coo_csr, 2*this_aero%tncoordinates, this%sparse_T2n%csr_values, &
                  this%sparse_T2n%csr_columns, this%sparse_T2n%csr_rowIndices, this%sparse_T2n%ncoords, & 
                  this%sparse_T2n%coo_values(1:this%sparse_T2n%ncoords), this%sparse_T2n%coo_rows(1:this%sparse_T2n%ncoords), this%sparse_T2n%coo_cols(1:this%sparse_T2n%ncoords+1), info)

    deallocate(this%sparse_Tn%coo_values)
    deallocate(this%sparse_Tn%coo_cols)
    deallocate(this%sparse_Tn%coo_rows)
    deallocate(this%sparse_T2n%coo_values)
    deallocate(this%sparse_T2n%coo_cols)
    deallocate(this%sparse_T2n%coo_rows)

    return
  end subroutine constructor
  
!< subroutine amplitude_material_loads
  subroutine amplitude_material_loads(this, this_structure, q1)
    implicit none

    class(model_fsi), intent(inout) :: this
    class(model_structure), intent(inout) :: this_structure
    
    real(kind = 8), allocatable, intent(in) :: q1(:)
    integer :: i
        
    this%ml_aero_t = 0.0d0
    
    do i = 1, this_structure%nnodes12
      !< determining the material load terms      
      call this_structure%nodes12(i)%materialterms(q1(this_structure%nodes12(i)%coordinates), this%fae(this_structure%nodes12(i)%coordinates))
      this%ml_aero_t(6*(i-1)+1:6*(i-1)+6) = this_structure%nodes12(i)%ml_t
    end do
    
    !< still to do for node6 !
    !do i = 1, this_structure%nnodes6
    !  !< determining the material load terms      
    !  call this_structure%nodes6(i)%materialterms(q1(this_structure%nodes6(i)%coordinates), this%fae(this_structure%nodes6(i)%coordinates))
    !  this%ml_aero_t(6*(i-1)+1:6*(i-1)+6) = this_structure%nodes6(i)%ml_t
    !end do
    
  end subroutine amplitude_material_loads
  
!< Subroutine to read input files
  subroutine readinginputs(this, this_structure, this_aero)
    implicit none
    
    class(model_fsi), intent(inout) :: this
    class(model_aero), intent(inout) :: this_aero
    class(model_structure), intent(inout) :: this_structure
    
    integer :: io_error
    integer :: i, j, k, i1
    integer, parameter :: linejump = 3
    integer :: isurface, innodes, nstruct
    character(100) :: strlfsort, strvinfsort, strfilename, stroldfilenameaero, stroldfilenamestructure, strstepsort, strfsitype, str_searchradius_filename, strfilename_windfield, str_fsi_type
    character(100) :: strfsiset
    real(kind = 8) :: searchradius
    integer :: icheck_static, icheck_dynamic, icheck_weak, icheck_strong, ifsiset, icheck_kinematic
    integer :: io

!< reading input file of vortex sheet
    print*, 'Parsing FSI-simulation input file'
    open(unit = 1000, file = 'simulationinput_fsi.txt', status = 'old', iostat = io_error)
      if (io_error == 0) then
        !< reading filenames
        do i = 1, linejump
          read(1000, *)
        end do
        read(1000, *) strfilename, stroldfilenameaero, stroldfilenamestructure
        if (allocated(this_aero%output_filename))      deallocate(this_aero%output_filename)
        if (allocated(this_structure%output_filename)) deallocate(this_structure%output_filename)
        allocate(this%output_filename, source           =  trim(adjustl(strfilename)))
        allocate(this_aero%output_filename, source      =  trim(adjustl(strfilename)))
        allocate(this_structure%output_filename, source =  trim(adjustl(strfilename)))
        
        if (trim(adjustl(stroldfilenameaero)) .ne. 'none' .and. trim(adjustl(stroldfilenameaero)) .ne. '!!') then
          this_aero%boolean_oldsim = .TRUE.
          allocate(this_aero%oldsimfilename, source =  trim(adjustl(stroldfilenameaero)))
        end if

        if (trim(adjustl(stroldfilenamestructure)) .ne. 'none' .and. trim(adjustl(stroldfilenamestructure)) .ne. '!!') then
          this_structure%boolean_oldsim = .TRUE.
          allocate(this_structure%oldsimfilename, source =  trim(adjustl(stroldfilenamestructure)))
        end if
        
        if (trim(adjustl(stroldfilenameaero)) .eq. '!!') backspace(1000)
        if (trim(adjustl(stroldfilenamestructure)) .eq. '!!') backspace(1000)
        
        !< reading aero simulation settings
        do i = 1, linejump
          read(1000, *)
        end do
        read(1000, *) this_aero%totalt, this_aero%deltat, this_aero%cutoff, this%timef, strlfsort
        this_aero%nsteps = IDNINT((this_aero%totalt-this_aero%time)/this_aero%deltat) + 1
        allocate(this%lfsort, source =  trim(adjustl(strlfsort)))
        
        !< reading wind settings
        do i = 1, linejump
          read(1000, *)
        end do
        read(1000, *) strvinfsort, this_aero%airflow%density, this_aero%airflow%intensity, this_aero%airflow%duration, this_aero%airflow%dir
        if (allocated(this_aero%airflow%sort)) deallocate(this_aero%airflow%sort)
        allocate(this_aero%airflow%sort, source =  trim(strvinfsort))
        
        !< reading structure simulation settings
        do i = 1, linejump
          read(1000, *)
        end do
        
        if (allocated(this_structure%settings)) deallocate(this_structure%settings)
        allocate(this_structure%settings(1))
        read(1000,  '(A)',  iostat = io_error) strfsiset

        icheck_kinematic = index(trim(adjustl(strfsiset)), 'kinematic')
        icheck_static  = index(trim(adjustl(strfsiset)), 'static')
        icheck_dynamic = index(trim(adjustl(strfsiset)), 'dynamic')
        icheck_weak    = index(trim(adjustl(strfsiset)), 'weak')
        icheck_strong  = index(trim(adjustl(strfsiset)), 'strong')
        
        if (icheck_kinematic .gt. 0) then
          allocate(this%str_fsi_type, source =  trim(adjustl(strfsiset)))
          allocate(this_structure%settings(1)%simutype, source = this%str_fsi_type)
        elseif (icheck_weak .gt. 0)  then
          allocate(this%str_fsi_type, source =  trim(adjustl(strfsiset(icheck_weak:len(strfsiset)))))
          ifsiset = icheck_weak-1
        elseif (icheck_strong .gt. 0) then
          allocate(this%str_fsi_type, source =  trim(adjustl(strfsiset(icheck_strong:len(strfsiset)))))
          ifsiset = icheck_strong-1
        else
          allocate(this%str_fsi_type, source = 'strong')
          ifsiset = len(strfsiset)
        end if

        if (icheck_static .gt. 0)   allocate(this_structure%settings(1)%simutype, source = trim(adjustl(strfsiset(icheck_static:ifsiset))))
        if (icheck_dynamic .gt. 0)  allocate(this_structure%settings(1)%simutype, source = trim(adjustl(strfsiset(icheck_dynamic:ifsiset))))
        
        if (this_structure%settings(1)%simutype .eq. 'dynamic') this_structure%settings(1)%i_simutype = 0
        if (this_structure%settings(1)%simutype .eq. 'static')  this_structure%settings(1)%i_simutype = 1

        read(1000, *) this_structure%settings(1)%deltat, this_structure%settings(1)%tolerance, this_structure%settings(1)%iterationlimit, this_structure%settings(1)%gravityflag
        read(1000, *, iostat=io) this_structure%settings(1)%output_flag, this%flag_dfae

        if (this%str_fsi_type .eq. 'strong') this%boolean_strongfsi = .TRUE.
        if (this%str_fsi_type .eq. 'weak') then
          this%boolean_strongfsi = .FALSE.
          this%flag_dfae = 0
        end if
        
        !< reading gravitiy vector
        do i = 1, linejump
          read(1000, *)
        end do
        read(1000, *) this_structure%gravity(:)
        
      !< Here external airflow data from file, if this_aero%airflow%sort = file
        if(adjustl(trim(this_aero%airflow%sort)) .eq. 'file') then 
          do i = 1, linejump
            read(1000, *)
          end do
        
          read(1000, *) strfilename_windfield
          allocate(this_aero%airflow%strfilename_windfield, source = adjustl(trim(strfilename_windfield)))
        
          do i = 1, linejump
            read(1000, *)
          end do
          read(1000, *) this_aero%airflow%center
        else
          !< checking for correct airflow sort - if not correct defined, then boolean_abort_airflow = .TRUE.
          this_aero%airflow%vinfinity = this_aero%airflow%dir * this_aero%evaluateamplitude(this_aero%time)
          if (this_aero%boolean_abort_airflow) then
            print*,'Warning: continueing with constant flow field...'
            this_aero%airflow%sort = 'constant'
            this_aero%boolean_abort_airflow = .FALSE.
          end if
        end if
        
      else
        print*,'Warning: no simulation file for DeSiO-FSI available...'
        this%boolean_abort = .TRUE.
        return
      end if
    close(unit = 1000)
    
    !< read airflow informations, if available
    if (allocated(this_aero%airflow%strfilename_windfield)) then
      this_aero%airflow%sort = 'constant'
      call this_aero%airflow%getairflowDatafromfile()
      if (this_aero%airflow%boolean_abort) then
        this%boolean_abort = .TRUE.
        return
      end if
    end if
    
    print*, 'Parsing FSI input file'
    open(unit = 1000, file = 'fsi_input.txt', status = 'old', iostat = io_error)
      if (io_error == 0) then
        !< reading fsi simulation settings
        do i = 1, linejump
          read(1000, *)
        end do
        read(1000, *) this%nfsi
        allocate(this%fsi(this%nfsi))
        
        !< loop over fsi
        do i = 1, this%nfsi
          
          do j = 1, linejump
            read(1000, *)
          end do          
          
          read(1000, *,  iostat = io_error) strfsitype, this%fsi(i)%fluid_surface, nstruct, searchradius
          !< checking if fsi_search_radius input is single value for all nodes or given by file for each aero-nodes
          if (io_error .ne. 0) then
            backspace(1000)
            read(1000, *,  iostat = io_error) strfsitype, this%fsi(i)%fluid_surface, nstruct, str_searchradius_filename
            if (io_error .ne. 0) then
              print*, 'error: in fsi_input search radius not defined...'
            end if
            allocate(this%fsi(i)%str_searchradius_filename, source = trim(str_searchradius_filename))        
          end if
          allocate(this%fsi(i)%structural_body(nstruct))
          read(1000, *) this%fsi(i)%structural_body

          if ( adjustl(trim(strfsitype)) .eq. 'beam') then
            this%fsi(i)%strtype = 'beam'
          elseif ( adjustl(trim(strfsitype)) .eq. 'shell') then
            this%fsi(i)%strtype = 'shell'
          elseif ( adjustl(trim(strfsitype)) .eq. 'rigid_body') then
            this%fsi(i)%strtype = 'rigid_body'
          end if
          
          if (allocated(this%fsi(i)%str_searchradius_filename)) then
            open(unit = 1001, file = this%fsi(i)%str_searchradius_filename, status = 'old', iostat = io_error)
              if (io_error .eq. 0 ) then
                do j = 1, linejump
                  read(1001, *)
                end do
                read(1001, *,  iostat = io_error) innodes
                do j = 1, linejump
                  read(1001, *)
                end do
                  
                allocate(this%fsi(i)%arr_searchradius(innodes,2))
                this%fsi(i)%arr_searchradius = 0.0d0
                !if (innodes .ne. this_aero%surfaces()) then
                do j = 1,innodes
                  this%fsi(i)%arr_searchradius(j,1) = j  
                  read(1001, *,  iostat = io_error) this%fsi(i)%arr_searchradius(j,2)
                end do
              else
                print*, 'error: in fsi_input search radius file ', this%fsi(i)%str_searchradius_filename
                this%boolean_abort = .TRUE.
                return
              end if
              
            close(unit = 1001)
          else
            innodes = this_aero%surfaces(this%fsi(i)%fluid_surface)%nnodes
            allocate(this%fsi(i)%arr_searchradius(innodes,2))
            this%fsi(i)%arr_searchradius = 0.0d0
            do j = 1,innodes
              this%fsi(i)%arr_searchradius(j,1) = j  
              this%fsi(i)%arr_searchradius(j,2) = searchradius
            end do
          end if
        end do
      else
        print*,'Warning: no fsi input file for DeSiO-FSI available...'
        this%boolean_abort = .TRUE.
        return
      end if
    close(unit = 1000)
    
    return
  end subroutine readinginputs

!< Subroutine to set and actualize kinematics of aerodynamical grid
  subroutine actualize(this, this_aero, q1, tncoordinates, v1_opt)
    implicit none

    class(model_aero), intent(inout) :: this_aero
    class(model_fsi),  intent(in)    :: this
    real(kind = 8),    intent(in)    :: q1(:)
    real(kind = 8),    intent(in), optional :: v1_opt(:)
    integer,           intent(in)    :: tncoordinates
    
    this_aero%qs_t = 0.0d0
    call mkl_dcsrmv('N', this_aero%tncoordinates, tncoordinates, 1.0d0, 'G  F  ', &
        this%sparse_Tn%csr_values, this%sparse_Tn%csr_columns, &
        this%sparse_Tn%csr_rowIndices(1:size(this%sparse_Tn%csr_rowIndices) - 1), &
        this%sparse_Tn%csr_rowIndices(2:size(this%sparse_Tn%csr_rowIndices)), &
        q1, 0.0d0, this_aero%qs_t)
    
        this_aero%qs_t = this_aero%qs_t + this%an
    
    if (present(v1_opt)) then
      this_aero%vs_t = 0.0d0
      call mkl_dcsrmv('N', this_aero%tncoordinates, tncoordinates, 1.0d0, 'G  F  ', &
          this%sparse_Tn%csr_values, this%sparse_Tn%csr_columns, &
          this%sparse_Tn%csr_rowIndices(1:size(this%sparse_Tn%csr_rowIndices) - 1), &
          this%sparse_Tn%csr_rowIndices(2:size(this%sparse_Tn%csr_rowIndices)), &
          v1_opt, 0.0d0, this_aero%vs_t)
    end if
    
    return
  end subroutine actualize
  
!< function for load factor, to scale the impulsive start situation
  function calc_load_factor(this, time) result(load_factor)
    implicit none 
    
    class(model_fsi) :: this
    real(kind = 8) :: time, load_factor
    
    load_factor = 0.0d0
    
    if (time .ge. this%timef) then
      load_factor = 1.0d0
    else
      if (this%lfsort .eq. 'cubic') then
        load_factor = (3.0d0*time**2)/this%timef**2 - (2.0d0*time**3)/this%timef**3
      elseif (this%lfsort .eq. 'tanh') then
        load_factor = 0.5d0*(exp(4.0d0*pi* (2.0d0*time-this%timef)/this%timef ) - 1.0d0 ) / (exp(4.0d0*pi* (2.0d0*time-this%timef)/this%timef ) + 1.0d0 ) - 0.5d0*(exp(-4.0d0*pi) - 1.0d0 ) / (exp(-4.0d0*pi) + 1.0d0 )
      elseif (this%lfsort .eq. 'linear') then
        load_factor = time/this%timef
      elseif (this%lfsort .eq. 'constant') then
        load_factor = 1.0d0
      elseif (this%lfsort .eq. 'null') then
        load_factor = 0.0d0
      end if
    endif

    return
  end function calc_load_factor
  
  end module class_model_fsi
