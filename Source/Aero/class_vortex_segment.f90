module class_vortex_segment

  use my_constants_aero, only: pi
  use my_types_aero, only: tangent_matrices3_3x3
  use my_math_aero, only: cross
  use class_node3_aero
  use my_vlm_functions, only: segment_velocity, tangent_segment_velocity
  
  implicit none
    
  type :: vortex_segment

     type(node3_aero) :: nodes(2)
     integer :: connectivity(2), connectivity_to_bounded_sheet(2) !< nodes
     integer :: adjacency(2) !< rings
     integer :: wakeedge = 0 !< if 1 the segment is an edge from which a wake is shed.
     integer :: first_row_wakeedge = 0
     integer :: surfaceedge = 0 !< if 1 the segment is an edge of the related surface.
     real :: time = 0.0d0 !< time step for diffusion model, the higher, the more segment circulation is reduced
     real(kind = 8) :: circulation = 0.0d0
     real(kind = 8) :: oldcirculation = 0.0d0
     real(kind = 8) :: cutoff = 0.0d0
     real(kind = 8) :: director(3)
     real(kind = 8) :: length
     real(kind = 8) :: sign_val = 0.0d0
     integer :: indicesq(6), globalID = 0
     
   contains
     
     procedure :: setgeometry => vortexsegmentsetgeometry
     procedure :: velocity => vortexsegmentvelocity
     procedure :: tangent_velocity => tangent_vortexsegmentvelocity
   end type vortex_segment

  contains

!< subroutine to set and actualize segment geometry
  subroutine vortexsegmentsetgeometry(this, q_t)
    implicit none
    
    class(vortex_segment), intent(inout) :: this
    real(kind = 8), intent(in) :: q_t(6)

    ! segment nodes
    this%nodes(1)%q_t(:) = q_t(1:3)
    this%nodes(2)%q_t(:) = q_t(4:6)

    ! tangent and length (not used)
    this%director = this%nodes(2)%q_t-this%nodes(1)%q_t
    this%length   = norm2(this%director)
    
    return
  end subroutine vortexsegmentsetgeometry
  
!< Function to calculate segment velocity  
  function vortexsegmentvelocity(this, targetpoint) result(v)
    implicit none
    
    class(vortex_segment) :: this
    real(kind = 8), intent(in) :: targetpoint(3)
    real(kind = 8) :: v(3)
    
    !< considering vortex-core growth model according to Lamb-Ossen (without viscous effects)
    v = segment_velocity(this%nodes(1)%q_t, this%nodes(2)%q_t, this%circulation, this%cutoff, targetpoint)
    
    return
  end function vortexsegmentvelocity
  
!< Tangent matrix due to linearized vortex-induced segment velocities of unbounded surface unit circulation
  function tangent_vortexsegmentvelocity(this, targetpoint) result(tangent_matrices)
    implicit none
    class(vortex_segment) :: this
    type(tangent_matrices3_3x3) :: tangent_matrices
    real(kind = 8), intent(in) :: targetpoint(3)
    real(kind = 8) :: kvseg(3,6)
    real(kind = 8) :: cutoff
    
    !< tangent matrix of vortex segment with unit circulation, depending to distances r1 and r2
    kvseg = tangent_segment_velocity(this%nodes(1)%q_t, this%nodes(2)%q_t, this%circulation, this%cutoff, targetpoint)
    
    !< reformulate tangent matrix for considering for ring nodes coordinates and nodal coordinates
    tangent_matrices%kvs1 = kvseg(1:3,1:3) + kvseg(1:3,4:6)
    tangent_matrices%kvs2 = -kvseg(1:3,1:3)
    tangent_matrices%kvs3 = -kvseg(1:3,4:6)
    
    return
  end function tangent_vortexsegmentvelocity
  
  end module class_vortex_segment
