module class_beam_element

  use my_constants_structure, only: i_33
  use my_math_structure, only: eye, outer
  use class_beam_element_dissipation
  implicit none
  
  type :: beam_element

     integer :: property
     integer :: connectivity(2)
     real(kind = 8), dimension(3, 3) :: cgg, ckk, cgk, cmass
     real(kind = 8) :: q_0(24)
     real(kind = 8) :: deltat
     integer :: simutype
     logical :: flag_kgeo_on, flag_kmat_on
     real(kind = 8) :: alpha_v, alpha_s
     real(kind = 8) :: m(24, 24)
     real(kind = 8) :: mhat(24,24)
     real(kind = 8) :: penergy
     real(kind = 8) :: denergy
     real(kind = 8) :: fqint(24)
     real(kind = 8) :: fqdyn(24)
     real(kind = 8) :: fv(24)
     real(kind = 8) :: kqq(24, 24), kqq_geo(24, 24), kqq_mat(24, 24)
     real(kind = 8) :: kvv(24, 24)
     real(kind = 8) :: kvq(24, 24)
     real(kind = 8) :: kqv(24, 24)
     real(kind = 8) :: stress_resultants(6)
     
     type(beam_element_dissipation) :: dissipation
   contains

     procedure :: internalterms => beaminternalterms
     procedure :: massmatrix => beammassmatrix
     
  end type beam_element
  
contains
  
  function  b1matrix(d1, d2, d3, phitick)
    
    implicit none
    
    real(kind = 8), dimension(3, 12) :: b1matrix
    real(kind = 8), dimension(3), intent(in) :: d1, d2, d3, phitick
    
    b1matrix(:, :) = 0.0d0
    
    b1matrix(1,  1: 3) = d1
    b1matrix(1,  4: 6) = phitick
    b1matrix(2,  1: 3) = d2
    b1matrix(2,  7: 9) = phitick
    b1matrix(3,  1: 3) = d3
    b1matrix(3, 10:12) = phitick
    
    return
    
  end function b1matrix
  
  function b2matrix(d1, d2, d3, d1tick, d2tick, d3tick)
    
    implicit none
    
    real(kind = 8), dimension(3, 18) :: b2matrix
    real(kind = 8), dimension(3), intent(in) :: d1, d2, d3, d1tick, d2tick, d3tick
    
    b2matrix(:, :) = 0.0d0
    
    b2matrix(1,  7: 9) =-d3tick
    b2matrix(1, 10:12) = d3
    b2matrix(1, 13:15) = d2tick
    b2matrix(1, 16:18) =-d2
    
    b2matrix(2,  1: 3) = d3tick
    b2matrix(2,  4: 6) =-d3
    b2matrix(2, 13:15) =-d1tick
    b2matrix(2, 16:18) = d1
    
    b2matrix(3,  1: 3) =-d2tick
    b2matrix(3,  4: 6) = d2
    b2matrix(3,  7: 9) = d1tick
    b2matrix(3, 10:12) =-d1

    b2matrix = b2matrix*0.5d0

    return

  end function b2matrix

  function q0matrix(n1, n2)

    implicit none

    real(kind = 8), dimension(12, 24) :: q0matrix
    real(kind = 8), intent(in) :: n1, n2

    q0matrix(:, :) = 0.0d0

    q0matrix( 1: 3,  1: 3) = n1*i_33
    q0matrix( 4: 6,  4: 6) = n1*i_33
    q0matrix( 7: 9,  7: 9) = n1*i_33
    q0matrix(10:12, 10:12) = n1*i_33

    q0matrix( 1: 3, 13:15) = n2*i_33
    q0matrix( 4: 6, 16:18) = n2*i_33
    q0matrix( 7: 9, 19:21) = n2*i_33
    q0matrix(10:12, 22:24) = n2*i_33

    return

  end function q0matrix
  
  function q1matrix(n1, n2, n1tick, n2tick)

    implicit none

    real(kind = 8), dimension(12, 24) :: q1matrix
    real(kind = 8), intent(in) :: n1, n2, n1tick, n2tick

    q1matrix(:, :) = 0.0d0

    q1matrix( 1: 3,  1: 3) = n1tick*i_33
    q1matrix( 4: 6,  4: 6) = n1*    i_33
    q1matrix( 7: 9,  7: 9) = n1*    i_33
    q1matrix(10:12, 10:12) = n1*    i_33

    q1matrix( 1: 3, 13:15) = n2tick*i_33
    q1matrix( 4: 6, 16:18) = n2    *i_33
    q1matrix( 7: 9, 19:21) = n2    *i_33
    q1matrix(10:12, 22:24) = n2    *i_33

    return

  end function q1matrix

  function q2matrix(n1, n2, n1tick, n2tick)

    implicit none

    real(kind = 8), dimension(18, 24) :: q2Matrix
    real(kind = 8), intent(in) :: n1, n2, n1Tick, n2Tick

    q2matrix(:, :) = 0.0d0

    q2matrix( 1: 3,  4: 6) = n1    *i_33
    q2matrix( 4: 6,  4: 6) = n1tick*i_33
    q2matrix( 7: 9,  7: 9) = n1    *i_33
    q2matrix(10:12,  7: 9) = n1tick*i_33
    q2matrix(13:15, 10:12) = n1    *i_33
    q2matrix(16:18, 10:12) = n1tick*i_33

    q2matrix( 1: 3, 16:18) = n2    *i_33
    q2matrix( 4: 6, 16:18) = n2tick*i_33
    q2matrix( 7: 9, 19:21) = n2    *i_33
    q2matrix(10:12, 19:21) = n2tick*i_33
    q2matrix(13:15, 22:24) = n2    *i_33
    q2matrix(16:18, 22:24) = n2tick*i_33
    
    return

  end function q2matrix

  function w11matrix(force)

    implicit none

    real(kind = 8), dimension(12, 12) :: w11matrix
    real(kind = 8), dimension(3), intent(in) :: force

    w11matrix(:, :) = 0.0d0

    w11matrix( 1: 3,  4: 6) = force(1)*i_33
    w11matrix( 1: 3,  7: 9) = force(2)*i_33
    w11matrix( 1: 3, 10:12) = force(3)*i_33
    w11matrix( 4: 6,  1: 3) = force(1)*i_33
    w11matrix( 7: 9,  1: 3) = force(2)*i_33
    w11matrix(10:12,  1: 3) = force(3)*i_33

    return

  end function w11matrix

  function w22matrix(moment)

    implicit none

    real(kind = 8), dimension(18, 18) :: w22matrix
    real(kind = 8), dimension(3), intent(in) :: moment

    w22Matrix(:, :) = 0.0d0

    w22matrix( 1: 3, 10:12) =-moment(3)*i_33
    w22matrix( 1: 3, 16:18) = moment(2)*i_33
    w22matrix( 4: 6,  7: 9) = moment(3)*i_33
    w22matrix( 4: 6, 13:15) =-moment(2)*i_33
    w22matrix( 7: 9,  4: 6) = moment(3)*i_33
    w22matrix( 7: 9, 16:18) =-moment(1)*i_33
    w22matrix(10:12,  1: 3) =-moment(3)*i_33
    w22matrix(10:12, 13:15) = moment(1)*i_33
    w22matrix(13:15,  4: 6) =-moment(2)*i_33
    w22matrix(13:15, 10:12) = moment(1)*i_33
    w22matrix(16:18,  1: 3) = moment(2)*i_33
    w22matrix(16:18,  7: 9) =-moment(1)*i_33

    w22matrix = w22matrix*0.5d0

    return

  end function w22matrix

  subroutine beaminternalterms(this, q_1, q_2_opt, v_1_opt, v_2_opt)
    implicit none
    class(beam_element) :: this
    real(kind = 8), dimension(24), intent(in) :: q_1
    real(kind = 8), dimension(24), optional,  intent(in) :: q_2_opt
    real(kind = 8), dimension(24), optional,  intent(in) :: v_1_opt
    real(kind = 8), dimension(24), optional,  intent(in) :: v_2_opt

    real(kind = 8), dimension(24) :: q_2, q_a
    real(kind = 8), dimension(3) :: deltaphi_0
    real(kind = 8) :: length_0
    real(kind = 8) :: n1, n2, n1tick, n2tick, jacobian, weight
    real(kind = 8), dimension(3) :: phi_0, d1_0, d2_0, d3_0, phitick_0, d1tick_0, d2tick_0, d3tick_0
    real(kind = 8), dimension(3) :: phi_1, d1_1, d2_1, d3_1, phitick_1, d1tick_1, d2tick_1, d3tick_1
    real(kind = 8), dimension(3) :: phi_2, d1_2, d2_2, d3_2, phitick_2, d1tick_2, d2tick_2, d3tick_2
    real(kind = 8), dimension(3) :: phi_a, d1_a, d2_a, d3_a, phitick_a, d1tick_a, d2tick_a, d3tick_a
    real(kind = 8), dimension(3) :: ggreek_1, kgreek_1, force_1, moment_1
    real(kind = 8), dimension(3) :: ggreek_2, kgreek_2, force_2, moment_2
    !
    real(kind = 8),  dimension(6,6)   :: c_beam             ! beam elasticity matrix
    real(kind = 8),  dimension(24,24) :: m                  ! consistent mass matrix
    real(kind = 8),  dimension(24,24) :: mhat               ! extended mass matrix == consist mass matrix and  identity matrix to consider v3 = d3dot
    
    real(kind = 8), dimension(6) :: strain_1, strain_2

    real(kind = 8), dimension(24) :: v_05                   ! velocity at n+0.5, sum of both

    !
    real(kind = 8), dimension(3) :: force_a, moment_a    
    real(kind = 8) :: b1_2(3, 12), b2_2(3, 18)
    real(kind = 8) :: b1_a(3, 12), b2_a(3, 18), c1_a(3, 24), c2_a(3, 24)
    real(kind = 8) :: q1(12, 24), q2(18, 24)
    real(kind = 8) :: m11(12, 12), m22(18, 18), m12(12, 18), m21(18, 12), w11(12, 12), w22(18, 18)
    
    real(kind = 8), dimension(24) :: v_1, v_2
    
    if (present(q_2_opt)) then
       q_2 = q_2_opt
    else
       q_2 = q_1
    end if

! Coordinates at qn+1/2    
    q_a = 0.5d0*(q_1+q_2)
    
    if (present(v_1_opt) .and. present(v_2_opt)) then
       v_1 = v_1_opt
       v_2 = v_2_opt
    else
       v_1(:) = 0.0d0
       v_2(:) = 0.0d0
    end if

    this%penergy   = 0.0d0
    this%denergy   = 0.0d0    
    this%fqint(:)  = 0.0d0
    this%fqdyn(:)  = 0.0d0
    this%fv(:)     = 0.0d0
    this%kqq(:, :) = 0.0d0
    this%kqq_geo(:, :) = 0.0d0
    this%kqq_mat(:, :) = 0.0d0
    this%kvv(:, :) = 0.0d0
    this%kqv(:, :) = 0.0d0
    this%kvq(:, :) = 0.0d0
    
! precalculation
    deltaphi_0 = this%q_0(13:15)-this%q_0( 1: 3)
    length_0 = dsqrt(dot_product(deltaphi_0, deltaphi_0))
        
    n1       = 0.5d0
    n2       = 0.5d0
    n1tick   =-1.0d0/length_0
    n2tick   = 1.0d0/length_0
    jacobian = 0.5d0*length_0
    weight   = 2.0d0
    
! Kinematic at t0
    phi_0     = n1    *this%q_0( 1: 3)+n2    *this%q_0(13:15)
    d1_0      = n1    *this%q_0( 4: 6)+n2    *this%q_0(16:18)
    d2_0      = n1    *this%q_0( 7: 9)+n2    *this%q_0(19:21)
    d3_0      = n1    *this%q_0(10:12)+n2    *this%q_0(22:24)

    phitick_0 = n1tick*this%q_0( 1: 3)+n2tick*this%q_0(13:15)
    d1tick_0  = n1tick*this%q_0( 4: 6)+n2tick*this%q_0(16:18)
    d2tick_0  = n1tick*this%q_0( 7: 9)+n2tick*this%q_0(19:21)
    d3tick_0  = n1tick*this%q_0(10:12)+n2tick*this%q_0(22:24)

! Kinematic at tn  
    phi_1     = n1    *q_1( 1: 3)+n2    *q_1(13:15)
    d1_1      = n1    *q_1( 4: 6)+n2    *q_1(16:18)
    d2_1      = n1    *q_1( 7: 9)+n2    *q_1(19:21)
    d3_1      = n1    *q_1(10:12)+n2    *q_1(22:24)

    phitick_1 = n1tick*q_1( 1: 3)+n2tick*q_1(13:15)
    d1tick_1  = n1tick*q_1( 4: 6)+n2tick*q_1(16:18)
    d2tick_1  = n1tick*q_1( 7: 9)+n2tick*q_1(19:21)
    d3tick_1  = n1tick*q_1(10:12)+n2tick*q_1(22:24)

! Kinematic at tn+1
    phi_2     = n1    *q_2( 1: 3)+n2    *q_2(13:15)
    d1_2      = n1    *q_2( 4: 6)+n2    *q_2(16:18)
    d2_2      = n1    *q_2( 7: 9)+n2    *q_2(19:21)
    d3_2      = n1    *q_2(10:12)+n2    *q_2(22:24)

    phitick_2 = n1tick*q_2( 1: 3)+n2tick*q_2(13:15)
    d1tick_2  = n1tick*q_2( 4: 6)+n2tick*q_2(16:18)
    d2tick_2  = n1tick*q_2( 7: 9)+n2tick*q_2(19:21)
    d3tick_2  = n1tick*q_2(10:12)+n2tick*q_2(22:24)

! Kinematic at Gau�point for coordinates at tn+1/2
    phi_a     = n1    *q_a( 1: 3)+n2    *q_a(13:15)
    d1_a      = n1    *q_a( 4: 6)+n2    *q_a(16:18)
    d2_a      = n1    *q_a( 7: 9)+n2    *q_a(19:21)
    d3_a      = n1    *q_a(10:12)+n2    *q_a(22:24)

    phitick_a = n1tick*q_a( 1: 3)+n2tick*q_a(13:15)
    d1tick_a  = n1tick*q_a( 4: 6)+n2tick*q_a(16:18)
    d2tick_a  = n1tick*q_a( 7: 9)+n2tick*q_a(19:21)
    d3tick_a  = n1tick*q_a(10:12)+n2tick*q_a(22:24)

! Beam-strain measures at tn
    ggreek_1(1) = dot_product(d1_1, phitick_1)-dot_product(d1_0, phitick_0) 
    ggreek_1(2) = dot_product(d2_1, phitick_1)-dot_product(d2_0, phitick_0)
    ggreek_1(3) = dot_product(d3_1, phitick_1)-dot_product(d3_0, phitick_0)    

    kgreek_1(1) = 0.5d0*((dot_product(d3_1, d2tick_1)-dot_product(d3_0, d2tick_0))-(dot_product(d2_1, d3tick_1)-dot_product(d2_0, d3tick_0)))
    kgreek_1(2) = 0.5d0*((dot_product(d1_1, d3tick_1)-dot_product(d1_0, d3tick_0))-(dot_product(d3_1, d1tick_1)-dot_product(d3_0, d1tick_0)))
    kgreek_1(3) = 0.5d0*((dot_product(d2_1, d1tick_1)-dot_product(d2_0, d1tick_0))-(dot_product(d1_1, d2tick_1)-dot_product(d1_0, d2tick_0)))

! Beam-strain measures at tn+1
    ggreek_2(1) = dot_product(d1_2, phitick_2)-dot_product(d1_0, phitick_0) 
    ggreek_2(2) = dot_product(d2_2, phitick_2)-dot_product(d2_0, phitick_0)
    ggreek_2(3) = dot_product(d3_2, phitick_2)-dot_product(d3_0, phitick_0)    

    kgreek_2(1) = 0.5d0*((dot_product(d3_2, d2tick_2)-dot_product(d3_0, d2tick_0))-(dot_product(d2_2, d3tick_2)-dot_product(d2_0, d3tick_0)))
    kgreek_2(2) = 0.5d0*((dot_product(d1_2, d3tick_2)-dot_product(d1_0, d3tick_0))-(dot_product(d3_2, d1tick_2)-dot_product(d3_0, d1tick_0)))
    kgreek_2(3) = 0.5d0*((dot_product(d2_2, d1tick_2)-dot_product(d2_0, d1tick_0))-(dot_product(d1_2, d2tick_2)-dot_product(d1_0, d2tick_0)))

! Stress resultants at (1) time n and (2) time n+1
    force_1  = matmul(this%cgg, ggreek_1)+matmul(          this%cgk , kgreek_1)
    moment_1 = matmul(this%ckk, kgreek_1)+matmul(transpose(this%cgk), ggreek_1) 

    force_2  = matmul(this%cgg, ggreek_2)+matmul(          this%cgk , kgreek_2)
    moment_2 = matmul(this%ckk, kgreek_2)+matmul(transpose(this%cgk), ggreek_2) 

! Elasticity matrix of the beam element    
    c_beam(1:3,1:3) = this%cgg
    c_beam(1:3,4:6) = this%cgk
    c_beam(4:6,1:3) = transpose(this%cgk)
    c_beam(4:6,4:6) = this%ckk
    
! Beam strains at time n and time n+1
    strain_1(1:3) = ggreek_1(1:3)
    strain_1(4:6) = kgreek_1(1:3)
    
    strain_2(1:3) = ggreek_2(1:3)
    strain_2(4:6) = kgreek_2(1:3)

! Mass matrices    
    m    = this%m           ! consistent mass matrix
    mhat = this%mhat        ! modified mass matrix to consider d3_dot = w3

! Velocity at tn+1/2
    v_05 = (v_1 + v_2)*0.5d0
    
! Initialization of dissipation
    this%dissipation%s_diss(:)      = 0.0d0
    this%dissipation%v_diss(:)      = 0.0d0
    this%dissipation%C_diss(:,:)    = 0.0d0
    this%dissipation%kvv_diss(:,:)  = 0.0d0
    
    if (this%simutype == 0 .or. this%simutype == 3) then
      if (this%alpha_s .ne. 0.0d0 .or. this%alpha_v .ne. 0.0d0) then
        call this%dissipation%evaluate(strain_1, strain_2, v_1, v_2, c_beam, m, this%alpha_s, this%alpha_v)
      end if
    end if 

! Calculating internal forces and stiffness matrices at tn+1/2 at Gausspoints
    force_a  = 0.5d0*(force_1  + force_2)  + this%dissipation%s_diss(1:3)
    moment_a = 0.5d0*(moment_1 + moment_2) + this%dissipation%s_diss(4:6)
    
    b1_2 = b1matrix(d1_2, d2_2, d3_2, phitick_2)                      ! B tn+1 Gamma
    b2_2 = b2matrix(d1_2, d2_2, d3_2, d1tick_2, d2tick_2, d3tick_2)   ! B tn+1 Kappa

    b1_a = b1matrix(d1_a, d2_a, d3_a, phitick_a)                      ! B-Operator tn+0.5 Gamma
    b2_a = b2matrix(d1_a, d2_a, d3_a, d1tick_a, d2tick_a, d3tick_a)   ! B-Operator tn+0.5 Kappa
    
    q1 = q1matrix(n1, n2, n1tick, n2tick)                             ! Ansatzfunctions for longitudinal strain 
    q2 = q2matrix(n1, n2, n1tick, n2tick)                             ! Ansatzfunctions for curvature

    c1_a = matmul(b1_a, q1)                                           ! Variation Gamma
    c2_a = matmul(b2_a, q2)                                           ! Variation Kappa

! B^T_{n+1/2}*c_beam*B_{n+1} for E_{n+1/2} - dissipation%C * 2.0d0 because of evaluation at E_{n+1}
    m11 = matmul(transpose(b1_a), matmul(this%cgg            + 2.0d0*this%dissipation%C_diss(1:3,1:3) , b1_2))
    m12 = matmul(transpose(b1_a), matmul(this%cgk            + 2.0d0*this%dissipation%C_diss(1:3,4:6) , b2_2))
    m21 = matmul(transpose(b2_a), matmul(transpose(this%cgk) + 2.0d0*this%dissipation%C_diss(4:6,1:3) , b1_2))
    m22 = matmul(transpose(b2_a), matmul(this%ckk            + 2.0d0*this%dissipation%C_diss(4:6,4:6) , b2_2))
    
    w11 = w11matrix(force_a)
    w22 = w22matrix(moment_a)
    
    this%kqq_mat = 0.0d0
    this%kqq_geo = 0.0d0
    
    if (this%flag_kgeo_on .eqv. .True.) then
      this%kqq_geo  = weight*jacobian*(matmul(transpose(q1), matmul(w11, q1)) + matmul(transpose(q2), matmul(w22, q2)))*0.5d0
    end if
    if (this%flag_kmat_on .eqv. .TRUE.) then
      this%kqq_mat  = weight*jacobian*(matmul(transpose(q1), matmul(m11, q1)) + matmul(transpose(q2), matmul(m22, q2)) + &
                                       matmul(transpose(q1), matmul(m12, q2)) + matmul(transpose(q2), matmul(m21, q1)))*0.5d0
    end if
    
    this%kqq    = this%kqq_geo + this%kqq_mat
    this%kqv    = m/this%deltat
    this%kvv    = -0.5d0*mhat - matmul(m,this%dissipation%kvv_diss)
    this%kvq    = mhat/this%deltat
    
    this%fqint  = weight*jacobian*(matmul(transpose(c1_a), force_a) + matmul(transpose(c2_a), moment_a)) ! = BQF = BN
    this%fqdyn  = matmul(m, (v_2-v_1)/this%deltat)
    this%fv     = matmul(mhat, (q_2-q_1)/this%deltat - v_05) - matmul(m, this%dissipation%v_diss)

    this%penergy = weight*jacobian*0.5d0*(dot_product(ggreek_2, force_2)+dot_product(kgreek_2, moment_2))
    
    this%stress_resultants(1:3) = force_2
    this%stress_resultants(4:6) = moment_2
    
    return
  end subroutine beaminternalterms

  function beammassmatrix(this) result(massmatrix)
   
    implicit none

    class(beam_element), intent(inout) :: this
    
    real(kind = 8), dimension(24, 24) :: massmatrix
    real(kind = 8), dimension(24, 24) :: mhat

    integer :: i
    real(kind = 8), dimension(3) :: deltaphi_0
    real(kind = 8) :: length_0
    real(kind = 8), dimension(2) :: n1, n2, weigth
    real(kind = 8) :: jacobian
    real(kind = 8), dimension(12, 24) :: q0
    real(kind = 8), dimension(12, 12) :: w00
    

    deltaphi_0 = this%q_0(13:15)-this%q_0( 1: 3)

    length_0 = dsqrt(dot_product(deltaphi_0, deltaphi_0))
    
    n1(1) = 0.788675134594813d0
    n1(2) = 0.211324865405187d0
    n2(1) = n1(2)
    n2(2) = n1(1)
    jacobian = 0.5d0*length_0
    weigth(:) = 1.0d0
       
    w00(:, :) = 0.0d0

    w00(1:3, 1:3) = this%cmass(1, 1)*i_33
    w00(1:3, 4:6) = this%cmass(1, 2)*i_33
    w00(1:3, 7:9) = this%cmass(1, 3)*i_33

    w00(4:6, 1:3) = this%cmass(2, 1)*i_33
    w00(4:6, 4:6) = this%cmass(2, 2)*i_33
    w00(4:6, 7:9) = this%cmass(2, 3)*i_33
    
    w00(7:9, 1:3) = this%cmass(3, 1)*i_33
    w00(7:9, 4:6) = this%cmass(3, 2)*i_33
    w00(7:9, 7:9) = this%cmass(3, 3)*i_33
  
    massmatrix(:, :) = 0.0d0

    this%m(:, :)   = 0.0d0
    this%mhat(:,:) = 0.0d0
    
    ! 2 point Gauss integration 
    do i = 1, 2 

       q0 = q0matrix(n1(i), n2(i))
       
       massmatrix = massmatrix+weigth(i)*matmul(transpose(q0), matmul(w00, q0))
       
    end do

    massmatrix =  massmatrix*jacobian
    mhat(:,:) = massmatrix(:,:)

! Due to momentum conservation of v3=d3_dot, the ones are included into the corresponding position of the consistent mass matrix
    mhat(10:12,10:12) = eye(3)      
    mhat(22:24,22:24) = eye(3)
    
    this%m = massmatrix
    this%mhat = mhat
    


    return
    
  end function beammassmatrix
  
end module class_beam_element
