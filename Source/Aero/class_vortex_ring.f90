module class_vortex_ring
  use my_constants_aero, only: pi, i_33
  use my_math_aero, only: cross, skew, outer
  use class_node3_aero
  use my_vlm_functions, only: segment_velocity, tangent_segment_velocity
  implicit none
  
  type :: vortex_ring
          
     type(node3_aero) :: nodes(4)
     integer :: connectivity(4) ! nodes
     integer :: adjacency(4) ! segments
     integer :: surfaceedge(4)
     integer :: wakeedge(4)
     
     real(kind = 8) :: circulation = 0.0d0, oldcirculation = 0.0d0
     real(kind = 8) :: cutoff = 0.0d0
     real(kind = 8) :: area
     real(kind = 8) :: center(3)
     real(kind = 8) :: normal(3), lvec(12)
     real(kind = 8) :: vmean(3)
     real(kind = 8) :: vsurf(3) = 0.0d0
     real(kind = 8) :: vsurfaces(3) = 0.0d0
     real(kind = 8) :: deltav(3)
     real(kind = 8) :: deltaps, deltaps_x, deltaps_w
     real(kind = 8) :: fae(3), fae_x(3), fae_w(3)
     real(kind = 8) :: vwakes(3) = 0.0d0
     real(kind = 8) :: Tcp(3,12) = 0.0d0                                                    ! transformation matrix for transformation of control point values to node values
     real(kind = 8) :: knormal(3,12) = 0.0d0, karea(12)     = 0.0d0, klvec(12,12) = 0.0d0
     real(kind = 8) :: kvsurf(3,12)  = 0.0d0, kvwakes(3,12) = 0.0d0
     real(kind = 8) :: vinfinity(3) = 0.0d0
     
     real(kind = 8), allocatable :: Tn(:,:)                                                 ! transformation matrix for transformation of aero node to structural node values
     real(kind = 8), allocatable :: avector(:)
     real(kind = 8), allocatable :: kqq_temp(:,:)
     integer, allocatable :: indicesq_st_global(:), indicesv_st_global(:)
     integer :: indicesv(12), indicesq(12), indicesp(3), indicesc
     
   contains
     procedure :: setgeometry => vortexringsetgeometry          
     procedure :: velocity => vortexringvelocity
     procedure :: setsurfvelocity => vortexsurfvelocity
     procedure :: aerodynamicloads => vortexringloads
     procedure :: getvmean
     procedure :: getdeltav
     procedure :: tangent_vortexringvelocity
     procedure :: tangent_vortexringgeometry
     procedure :: tangent_deltav
  end type vortex_ring

  contains

!< Subroutine for calculating aerodynamic forces  
  subroutine vortexringloads(this, circulations, segmentscirculations, oldsegmentscirculations, tnrings, deltat, tcounter, density, boolean_unsteady_term)
    implicit none 
    class(vortex_ring), intent(inout) :: this
    
    real(kind = 8), intent(in) :: circulations(:), segmentscirculations(:), oldsegmentscirculations(:), deltat, density
    integer, intent(in) :: tnrings, tcounter 
    logical, intent(in) :: boolean_unsteady_term
    real(kind = 8) :: circulation_dot, circulation_dot_x, circulation_dot_w
    
    this%vmean        = this%getvmean(circulations, tnrings)
    this%deltav       = this%getdeltav(segmentscirculations, oldsegmentscirculations)
    
    circulation_dot    = 0.0d0
    if (boolean_unsteady_term) then
        if (tcounter .gt. 0) then
          circulation_dot   = ( this%circulation - this%oldcirculation ) / (deltat)
        end if
    end if
    
    !< complete aerodynamic force
    this%deltaps = density * ( dot_product(this%vmean - this%vsurf, this%deltav) - circulation_dot )
    this%fae     = this%deltaps * this%normal * this%area

    return
  end subroutine vortexringloads
  
!< Function to calculate vmean of fluid
  function getvmean(this, circulations, tnrings) result(vmean)
    implicit none 
    class(vortex_ring), intent(inout) :: this
    
    real(kind = 8), intent(in ):: circulations(:)
    integer, intent(in):: tnrings
    real(kind = 8) :: vsurfaces(3), vmean(3)
    integer :: k

    !< velocity due to unit surface circulations already calculated in evaluateamatrix: this%surfaces(i)%rings(j)%avector
    vsurfaces(:) = 0.0d0
    do k = 1,tnrings
      vsurfaces = vsurfaces + this%avector(3*(k-1)+1:3*(k-1)+3)*circulations(k)
    end do
    this%vsurfaces = vsurfaces
    
    !< velocity due to wake segments already calculated in evaluaterhs: this%surfaces(i)%rings(j)%vwakes
    !< absolute velocity at ring control point: only velocity of the fluid -> see derivation based on Bernoulli equation
    vmean = this%vinfinity + this%vsurfaces + this%vwakes
    
    return
  end function getvmean
  
!< Function to calculate deltav
  function getdeltav(this, segmentcirculations, oldsegmentscirculations) result(deltav)
    implicit none 
    class(vortex_ring), intent(inout) :: this
    real(kind = 8), intent(in):: segmentcirculations(:), oldsegmentscirculations(:)
    real(kind = 8) :: deltav(3), circvec(3), circ
    real(kind = 8) :: factor
    integer :: k
    
    !< determine vector of circulations of a ring for calculation of steady forces
    circvec = 0.0d0 
    do k = 1, 4
      ! differentiate between outer surface edge and inner surface edge. Segment circulation in a ring = of the global segment cirulation assigned in class_vortex_sheet (line 50)
        factor = 1.0d0
        if (this%surfaceedge(k) == 1) factor = 2.0d0
        if (this%wakeedge(k) == 1) then
          factor = 1.0d0
          circ = segmentcirculations(this%adjacency(k)) - oldsegmentscirculations(this%adjacency(k))
        else
          circ = factor*segmentcirculations(this%adjacency(k))
        end if
        circvec = circvec + 0.5d0*circ*this%lvec(3*(k-1)+1:3*(k-1)+3)
    end do
    
    !< deltav = -n x circvec / A
    deltav = - cross(this%normal, circvec)/this%area

    return
  end function getdeltav
  
!< Subroutine to set vortex ring geometry  
  subroutine vortexringsetgeometry(this, q_t)
    implicit none
    class(vortex_ring), intent(inout) :: this
    real(kind = 8), intent(in) :: q_t(12)
    real(kind = 8), dimension(3) :: e1, e2, d12, d21, d14, d13, d24, d41, d43, d23, L1, L2, h1, h2
    
    this%nodes(1)%q_t(:) = q_t( 1: 3)
    this%nodes(2)%q_t(:) = q_t( 4: 6)
    this%nodes(3)%q_t(:) = q_t( 7: 9)
    this%nodes(4)%q_t(:) = q_t(10:12)

    !< calculating center of ring surface
    this%center(:) = 0.0d0
    this%center = this%center + 0.25d0* ( this%nodes(1)%q_t + this%nodes(2)%q_t + this%nodes(3)%q_t + this%nodes(4)%q_t )
    
    !< calculating diagonal of ring surface
    d12 = this%nodes(1)%q_t-this%nodes(2)%q_t
    d13 = this%nodes(1)%q_t-this%nodes(3)%q_t
    d14 = this%nodes(1)%q_t-this%nodes(4)%q_t
    d21 = this%nodes(2)%q_t-this%nodes(1)%q_t
    d23 = this%nodes(2)%q_t-this%nodes(3)%q_t
    d24 = this%nodes(2)%q_t-this%nodes(4)%q_t
    d41 = this%nodes(4)%q_t-this%nodes(1)%q_t
    d43 = this%nodes(4)%q_t-this%nodes(3)%q_t
    
    !< calculating unit vectors of (trapezodial) ring surface
    h1 = d24-d13
    h2 = -(d13+d24)
    e1 = h1/norm2(h1)
    e2 = h2/norm2(h2)
    
    !< calculating plane normal of ring surface
    this%normal = cross(e1,e2)
    
    ! calculating lvec with respect to global segment orientation
    this%lvec(1:3)   = d43 ! segment 1
    this%lvec(4:6)   = d14 ! segment 2
    this%lvec(7:9)   = d12 ! segment 3
    this%lvec(10:12) = d23 ! segment 4

    !< calculating plane area using formula for triangular surfaces
    L1  = cross(d21,d41)
    L2  = cross(d43,d23)
    this%area = 0.5d0 * (norm2(L1) + norm2(L2))

    return
  end subroutine vortexringsetgeometry

!< Subroutine to calculate vortex ring velocity
  function vortexringvelocity(this, targetpoint) result(v)
    implicit none
    class(vortex_ring) :: this
    real(kind = 8), intent(in) :: targetpoint(3)
    real(kind = 8) :: v(3)
    integer :: i, j
    
    v(:) = 0.0d0
    do i = 1, 4
       j = i+1
       if (i == 4) then
        j = 1 
       end if
       v = v+segment_velocity(this%nodes(i)%q_t, this%nodes(j)%q_t, this%circulation, this%cutoff, targetpoint)
    end do
    
    return
  end function vortexringvelocity

!< Subroutine to set prescribed vortex ring velocity, i.e. of bounded vortex sheet
  subroutine vortexsurfvelocity(this, v_t)
    implicit none
    class(vortex_ring), intent(inout) :: this
    real(kind = 8), intent(in) :: v_t(12)
        
    this%nodes(1)%v_t(:) = v_t( 1: 3)
    this%nodes(2)%v_t(:) = v_t( 4: 6)
    this%nodes(3)%v_t(:) = v_t( 7: 9)
    this%nodes(4)%v_t(:) = v_t(10:12)

    this%vsurf(:) = 0.0d0
    
    !< calculating center of ring surface
    this%vsurf = 0.25d0*( this%nodes(1)%v_t + this%nodes(2)%v_t + this%nodes(3)%v_t + this%nodes(4)%v_t )
        
    return
  end subroutine vortexsurfvelocity
 
! =========================================================================================================
! TANGENT MATRICES
! =========================================================================================================
!< Tangent matrix of deltav
  subroutine tangent_deltav(this, segmentcirculations, oldsegmentscirculations, kdv1, indicessegmentcirc, factor_tvec)
    implicit none 
    class(vortex_ring), intent(inout) :: this
    real(kind = 8), intent(in):: segmentcirculations(:), oldsegmentscirculations(:)
    real(kind = 8), intent(inout) :: kdv1(3,12), factor_tvec(12)
    integer, intent(inout) :: indicessegmentcirc(4)
    real(kind = 8) :: circvec(3), circ
    real(kind = 8) :: kcircvec1(3,12)
    real(kind = 8) :: factor
    integer :: k
    
    !< determine vector of circulations of a ring for calculation of steady forces
    circvec = 0.0d0
    kcircvec1 = 0.0d0
    do k = 1, 4
      factor = 1.0d0
      if (this%surfaceedge(k) == 1) factor = 2.0d0
      if (this%wakeedge(k) == 1) then
        factor = 1.0d0
        circ = segmentcirculations(this%adjacency(k)) - oldsegmentscirculations(this%adjacency(k))
      else
        circ = factor*segmentcirculations(this%adjacency(k))
      end if
      circvec   = circvec   + 0.5d0*circ*this%lvec(3*(k-1)+1:3*(k-1)+3)
      kcircvec1 = kcircvec1 + 0.5d0*circ*this%klvec(3*(k-1)+1:3*(k-1)+3,1:12)
      factor_tvec(3*(k-1)+1:3*(k-1)+3) = 0.5d0*this%lvec(3*(k-1)+1:3*(k-1)+3)*factor
      indicessegmentcirc(k) = this%adjacency(k)
    end do
    
    !< D deltav = D (-n x circvec / A)
    kdv1 = matmul(skew(circvec),this%knormal)/this%area + outer(cross(this%normal,circvec),this%karea)/(this%area**2) - matmul(skew(this%normal),kcircvec1)/this%area 
    
    return
  end subroutine tangent_deltav
  
!< Tangent matrix due to linearized vortex-induced ring velocity of bounded surface unit circulation
  subroutine tangent_vortexringvelocity(this, kvs1, kvs2, targetpoint)
    implicit none
    class(vortex_ring) :: this
    real(kind = 8), intent(in) :: targetpoint(3)
    real(kind = 8), intent(inout) :: kvs1(3,3), kvs2(3,12)
    real(kind = 8) :: kvseg1(3,6), kvseg2(3,6), kvseg3(3,6), kvseg4(3,6)
    
    !< tangent matrix of vortex segment with unit circulation, depending to distances r1 and r2
    kvseg1 = tangent_segment_velocity(this%nodes(1)%q_t, this%nodes(2)%q_t, 1.0d0, this%cutoff, targetpoint)
    kvseg2 = tangent_segment_velocity(this%nodes(2)%q_t, this%nodes(3)%q_t, 1.0d0, this%cutoff, targetpoint)
    kvseg3 = tangent_segment_velocity(this%nodes(3)%q_t, this%nodes(4)%q_t, 1.0d0, this%cutoff, targetpoint)
    kvseg4 = tangent_segment_velocity(this%nodes(4)%q_t, this%nodes(1)%q_t, 1.0d0, this%cutoff, targetpoint)
    
    !< tangent matrix for control point of vortex ring
    kvs1 = 0.0d0 
    kvs1 = kvseg1(1:3,1:3) + kvseg1(1:3,4:6) + kvseg2(1:3,1:3) + kvseg2(1:3,4:6) + kvseg3(1:3,1:3) + kvseg3(1:3,4:6) + kvseg4(1:3,1:3) + kvseg4(1:3,4:6)
    
    !<  tangent depending to ring nodes coordinates
    kvs2 = 0.0d0
    kvs2(1:3,1:3)   = -(kvseg1(1:3,1:3) + kvseg4(1:3,4:6))
    kvs2(1:3,4:6)   = -(kvseg1(1:3,4:6) + kvseg2(1:3,1:3))
    kvs2(1:3,7:9)   = -(kvseg2(1:3,4:6) + kvseg3(1:3,1:3))
    kvs2(1:3,10:12) = -(kvseg3(1:3,4:6) + kvseg4(1:3,1:3))
     
    return
  end subroutine tangent_vortexringvelocity
 
!< Tangent matrix due to linearized area and area normal  
  subroutine tangent_vortexringgeometry(this)
    implicit none
    class(vortex_ring), intent(inout) :: this
    real(kind = 8), dimension(3) :: e1, e2, d12, d21, d14, d13, d24, d41, d43, d23, L1, L2, h1, h2
    real(kind = 8), dimension(3,12) :: B21, B24, B13, B41, B23, B43, B12, B14
    real(kind = 8) :: De1(3,12), De2(3,12)
    
    B21=0.0d0
    B24=0.0d0
    B13=0.0d0
    B41=0.0d0
    B23=0.0d0
    B43=0.0d0
    B12=0.0d0
    B14=0.0d0

    B12(1:3,1:3)   =  i_33
    B12(1:3,4:6)   = -i_33
    B13(1:3,1:3)   =  i_33
    B13(1:3,7:9)   = -i_33
    B14(1:3,1:3)   =  i_33
    B14(1:3,10:12) = -i_33
    B21(1:3,1:3)   = -i_33
    B21(1:3,4:6)   =  i_33
    B23(1:3,4:6)   =  i_33
    B23(1:3,7:9)   = -i_33
    B24(1:3,4:6)   =  i_33
    B24(1:3,10:12) = -i_33
    B41(1:3,1:3)   = -i_33
    B41(1:3,10:12) =  i_33
    B43(1:3,7:9)   = -i_33
    B43(1:3,10:12) =  i_33
    
    !< calculating diagonal of ring surface
    d12 = this%nodes(1)%q_t-this%nodes(2)%q_t
    d13 = this%nodes(1)%q_t-this%nodes(3)%q_t
    d14 = this%nodes(1)%q_t-this%nodes(4)%q_t
    d21 = this%nodes(2)%q_t-this%nodes(1)%q_t
    d23 = this%nodes(2)%q_t-this%nodes(3)%q_t
    d24 = this%nodes(2)%q_t-this%nodes(4)%q_t
    d41 = this%nodes(4)%q_t-this%nodes(1)%q_t
    d43 = this%nodes(4)%q_t-this%nodes(3)%q_t
    
    !< calculating unit vectors of (trapezodial) ring surface
    h1 = d24-d13
    h2 = -(d13+d24)
    e1 = h1/norm2(h1)
    e2 = h2/norm2(h2)
    
    !< calculating plane area using formula for triangular surfaces
    L1  = cross(d21,d41)
    L2  = cross(d43,d23)

    De1 = matmul(  i_33*norm2(h1) - outer(h1,h1)/norm2(h1),  (B24 - B13) ) / (dot_product(h1,h1))
    De2 = matmul(  i_33*norm2(h2) - outer(h2,h2)/norm2(h2), -(B24 + B13) ) / (dot_product(h2,h2))

    !< calculating tangent matrices for area normal, area and lvec
    this%knormal(1:3,1:12) = matmul(skew(e1),De2) - matmul(skew(e2),De1)
    this%karea(1:12)       = 0.50d0* ( matmul(L1,matmul(skew(d21),B41) - matmul(skew(d41),B21))/norm2(L1) + matmul(L2,matmul(skew(d43),B23) - matmul(skew(d23),B43))/norm2(L2) )
    this%klvec = 0.0d0
    !this%klvec(1:3,1:12)   = B12
    !this%klvec(4:6,1:12)   = B23
    !this%klvec(7:9,1:12)   = B43
    !this%klvec(10:12,1:12) = B14
    this%klvec(1:3,1:12)   = B43
    this%klvec(4:6,1:12)   = B14
    this%klvec(7:9,1:12)   = B12
    this%klvec(10:12,1:12) = B23

    return
  end subroutine tangent_vortexringgeometry
  
end module class_vortex_ring
