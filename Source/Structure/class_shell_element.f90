module class_shell_element

  use my_math_structure, only: cross, eye, diag, curva2curvb3matrix, curva2curvb6matrix
  use my_constants_structure, only: i_33
  use class_shell_element_dissipation
  use class_multi_layer_material

  implicit none
  
  type :: shell_element

     integer :: connectivity(4)
     integer :: coordinates(8) ! internal degrees of freedom due to the enhancements.
     real(kind = 8) :: q_0(32)
     real(kind = 8) :: deltat
     integer :: simutype
     logical :: flag_kgeo_on, flag_kmat_on
     real(kind = 8) :: e1_0(3)
     real(kind = 8) :: e2_0(3)
     real(kind = 8) :: e3_0(3)
     real(kind = 8) :: gmsup1_0(3)
     real(kind = 8) :: gmsup2_0(3)
     real(kind = 8) :: gmsup3_0(3)
     real(kind = 8) :: gmsub1_0(3)
     real(kind = 8) :: gmsub2_0(3)
     real(kind = 8) :: gmsub3_0(3)
     real(kind = 8) :: sqrtgm_0
     real(kind = 8) :: length_0
     real(kind = 8) :: penergy
     real(kind = 8) :: denergy
     real(kind = 8) :: fqint(32)
     real(kind = 8) :: fqdyn(32)
     real(kind = 8) :: fv(32)
     real(kind = 8) :: kqq(32, 32)
     real(kind = 8) :: kqq_geo(32, 32)
     real(kind = 8) :: kqq_mat(32, 32)
     
     real(kind = 8) :: kvv(32, 32)
     real(kind = 8) :: kvq(32, 32)
     real(kind = 8) :: kqv(32, 32)
     real(kind = 8) :: mass_matrix(24, 24)
     
     type(shell_element_dissipation) :: dissipation
     type(multi_layer_material) :: laminate
     
   contains

     procedure :: internalterms => shellinternalterms
     procedure :: massmatrix => shellmassmatrix
     procedure :: precalculation => shellprecalculation
     
  end type shell_element
  
contains

  subroutine shellprecalculation(this)

    implicit none
    
    class(shell_element), intent(inout) :: this

    real(kind = 8) :: d3zeta
    real(kind = 8) :: xi(3)
    real(kind = 8), dimension(4) :: n, d1n, d2n
    
    real(kind = 8), dimension(3, 24) :: nphi, ndir, d1nphi, d1ndir, d2nphi, d2ndir
    
    real(kind = 8), dimension(3) :: gmsub1_0, gmsub2_0, gmsub3_0
    
    real(kind = 8), dimension(3) :: ediag1_0, ediag2_0

    d3zeta = 0.5d0*this%laminate%thickness
    
    xi = [0.0d0, 0.0d0, 0.0d0]
    
    call unitsquareshapefunctions(n, d1n, d2n, xi)
    
    nphi   = nphimatrix(n)
    d1nphi = nphimatrix(d1n)
    d2nphi = nphimatrix(d2n)
    
    ndir   = ndirmatrix(n)
    d1ndir = ndirmatrix(d1n)
    d2ndir = ndirmatrix(d2n)             
    
    gmsub1_0 = matmul(d1nphi, this%q_0(1:24))
    
    gmsub2_0 = matmul(d2nphi, this%q_0(1:24))
    
    gmsub3_0 = d3zeta*matmul(ndir, this%q_0(1:24))
    
    this%sqrtgm_0 = dot_product(cross(gmsub1_0, gmsub2_0), gmsub3_0)
    
    this%gmsup1_0 = cross(gmsub2_0, gmsub3_0)/this%sqrtgm_0
    
    this%gmsup2_0 = cross(gmsub3_0, gmsub1_0)/this%sqrtgm_0
    
    this%gmsup3_0 = cross(gmsub1_0, gmsub2_0)/this%sqrtgm_0

    this%gmsub1_0 = gmsub1_0

    this%gmsub2_0 = gmsub2_0

    this%gmsub3_0 = gmsub3_0
        
    ediag1_0 = (this%q_0( 1: 3)-this%q_0(13:15))/norm2(this%q_0( 1: 3)-this%q_0(13:15)) 
    
    ediag2_0 = (this%q_0( 7: 9)-this%q_0(19:21))/norm2(this%q_0( 7: 9)-this%q_0(19:21))
    
    this%e1_0 = (ediag1_0-ediag2_0)/norm2(ediag1_0-ediag2_0)
    
    this%e2_0 = (ediag1_0+ediag2_0)/norm2(ediag1_0+ediag2_0)

    this%e3_0 = cross(this%e1_0, this%e2_0)
    
    return
    
  end subroutine shellprecalculation
  
  subroutine shellinternalterms(this, q_1, q_2_opt, v_1_opt, v_2_opt)
    
    implicit none

    class(shell_element) :: this
    real(kind = 8), intent(in) :: q_1(32)
    real(kind = 8), optional, intent(in) :: q_2_opt(32)
    real(kind = 8), optional, intent(in) :: v_1_opt(32)
    real(kind = 8), optional, intent(in) :: v_2_opt(32)
        
    real(kind = 8), dimension(32) :: q_2, q_a
    real(kind = 8), dimension(6) :: strain_1, strainc_1, straine_1, stress_1
    real(kind = 8), dimension(6) :: strain_2, strainc_2, straine_2, stress_2
    real(kind = 8), dimension(6) :: stress_a
    real(kind = 8), dimension(3, 3) :: curv2cart3
    real(kind = 8), dimension(6, 6) :: cscriptcurv, curv2cart6

    ! Number of Gaussian Points
    integer :: ngaussianpoints = 2  ! in plane integration
    integer :: ngaussianpointst = 2 ! thickness integration nGP = 1,2,3,4 or 5
     
    real(kind = 8), dimension(2) :: xigauss, wigauss
    real(kind = 8), dimension(5) :: xigausst, wigausst
    
    real(kind = 8) :: xi(3)
    real(kind = 8), dimension(4) :: n, d1n, d2n
    
    real(kind = 8) :: sqrtg_0
    real(kind = 8), dimension(3) :: gsub1_0, gsub2_0, gsub3_0
    real(kind = 8), dimension(3) :: gsup1_0, gsup2_0, gsup3_0
    
    real(kind = 8) :: zeta, d3zeta
        
    real(kind = 8), dimension(3, 24) :: nphi, ndir, d1nphi, d1ndir, d2nphi, d2ndir
    real(kind = 8), dimension(24, 24) :: m11, m22, m33, m23, m13, m12

    real(kind = 8) :: bbar_2(6, 24), bbartrans_2(24, 6)
    real(kind = 8) :: bbar_a(6, 24), bbartrans_a(24, 6), btotal_a(6, 32), btotal_2(6, 32)
    real(kind = 8) :: btotal_a_t(32, 6) ! transpose of btotal_a
    real(kind = 8) :: bhat(6, 8)
    real(kind = 8) :: w1(32, 32), w2(32, 32), kqq_mat(32,32), kqq_geo(32,32)
    
    integer :: i, j, k, l, info

    real(kind = 8) :: wijk_times_sqrtg_0
    
    real(kind = 8), dimension(32) :: v_1, v_2
    real(kind = 8), dimension(3, 32) :: nq
    real(kind = 8), dimension(32, 3) :: nqt ! transpose of nq
    real(kind = 8), dimension(32, 32) :: nq_square_dt ! (nq^T * nq) / dt
    
    real(kind = 8), dimension(3) :: qh_1, qh_2, vh_1, vh_2, vh_a
    real(kind = 8), dimension(6) :: deltastrain
    
    real(kind = 8) :: kee(8,8), RHS(8,24), kqq_sc(24,24), kqq_temp(32,32)           ! matrices for static condensation (simutype == 2, 3)
    real(kind = 8), dimension(:), allocatable :: ipiv
    
    !if nGP == 3; 
    if (ngaussianpointst ==1) then
        xigausst = 0.0d0*[0.0d0, 0.0d0,0.0d0,0.0d0,0.0d0]
        wigausst = 2.0d0*              [ 1.0d0, 0.0d0,0.0d0,0.0d0,0.0d0]
    elseif (ngaussianpointst ==2) then
        xigausst = 0.577350269189626d0*[-1.0d0, 1.0d0,0.0d0,0.0d0,0.0d0]
        wigausst =                     [ 1.0d0, 1.0d0,0.0d0,0.0d0,0.0d0]
    else if (ngaussianpointst ==3) then
        xigausst = [-0.77459667d0,0.0d0     ,0.77459667d0,0.0d0,0.0d0]
        wigausst = [ 0.555556d0  ,0.888889d0,0.555556d0  ,0.0d0,0.0d0]
    else if (ngaussianpointst ==4) then
        xigausst = [-0.8611363d0,-0.3399810436d0,0.3399810436d0,0.8611363d0,0.0d0]
        wigausst = [0.3478548d0,0.6521451549d0,0.6521451549d0,0.3478548d0,0.0d0]
    else if (ngaussianpointst ==5) then
        xigausst = [-0.9061798459d0,-0.5384693101d0,0d0,0.5384693101d0,0.9061798459d0]
        wigausst = [0.2369268851d0,0.4786286705d0,0.5688888889d0,0.4786286705d0,0.2369268851d0]
    end if
    
    xigauss = 0.577350269189626d0*[-1.0d0, 1.0d0]
    wigauss =                     [ 1.0d0, 1.0d0]

    d3zeta = this%laminate%thickness*0.5d0
    
    if (present(q_2_opt)) then

       q_2 = q_2_opt
       
    else

       q_2 = q_1
       
    end if

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
    this%kqq_mat(:, :) = 0.0d0
    this%kqq_geo(:, :) = 0.0d0
    this%kvv(:, :) = 0.0d0
    this%kqv(:, :) = 0.0d0
    this%kvq(:, :) = 0.0d0
    
    kqq_mat = 0.0d0
    kqq_geo = 0.0d0
    
    bbar_2(:, :) = 0.0d0
    bbar_a(:, :) = 0.0d0
    
    bhat(:, :) = 0.0d0

    btotal_a(:, :) = 0.0d0
    btotal_2(:, :) = 0.0d0
    
    w1(:, :) = 0.0d0
    w2(:, :) = 0.0d0
   
    nq(:, :) = 0.0d0
    
    do l = 1, this%laminate%nplys
       
       do k = 1, ngaussianpointst
       
            xi(3) = xigausst(k)
           
            zeta = this%laminate%boundaries(l)*(1.0d0-xi(3))*0.5d0+this%laminate%boundaries(l+1)*(1.0d0+xi(3))*0.5d0    
            
          ! Tranformation in isoparamterischen Raum (theta3)
            xi(3) = 2.0d0*zeta/this%laminate%thickness            
            
          do j = 1, ngaussianpoints

             do i = 1, ngaussianpoints ! based on isopararametric space [-1,1] for each layer
                 
                xi(1:2) = [xigauss(i), xigauss(j)]                                                                
                
                call unitsquareshapefunctions(n, d1n, d2n, xi)

                nphi   = nphimatrix(n)
                d1nphi = nphimatrix(d1n)
                d2nphi = nphimatrix(d2n)
             
                ndir   = ndirmatrix(n)
                d1ndir = ndirmatrix(d1n)
                d2ndir = ndirmatrix(d2n)

                nq(:, 1:24) = nphi+zeta*ndir

                qh_1 = matmul(nq, q_1)
                qh_2 = matmul(nq, q_2)
             
                vh_1 = matmul(nq, v_1)
                vh_2 = matmul(nq, v_2)             

                vh_a = 0.5d0*(vh_1+vh_2)
             
                gsub1_0 = matmul(d1nphi+zeta*d1ndir, this%q_0(1:24))

                gsub2_0 = matmul(d2nphi+zeta*d2ndir, this%q_0(1:24))

                gsub3_0 = d3zeta*matmul(ndir, this%q_0(1:24))
       
                sqrtg_0 = dot_product(cross(gsub1_0, gsub2_0), gsub3_0)

                gsup1_0 = cross(gsub2_0, gsub3_0)/sqrtg_0

                gsup2_0 = cross(gsub3_0, gsub1_0)/sqrtg_0

                gsup3_0 = cross(gsub1_0, gsub2_0)/sqrtg_0

				! Mapping between two isoparametric spaces / scaling of weights using the thickness ration
                wijk_times_sqrtg_0 = wigauss(i)*wigauss(j)*wigausst(k)*sqrtg_0*this%laminate%plys(l)%thickness/this%laminate%thickness
                
!TODO this call takes a lot of time. it seems that a lot of terms are constant at every time step and can be cached
                call mmatrices(m11, m22, m33, m23, m13, m12, xi, zeta, d3zeta, nphi, ndir, d1nphi, d1ndir, d2nphi, d2ndir)
             
                bbartrans_2(:, 1) = matmul(m11               , q_2(1:24))
                bbartrans_2(:, 2) = matmul(m22               , q_2(1:24))
                bbartrans_2(:, 3) = matmul(m33               , q_2(1:24))
                bbartrans_2(:, 4) = matmul(m23+transpose(m23), q_2(1:24))
                bbartrans_2(:, 5) = matmul(m13+transpose(m13), q_2(1:24))
                bbartrans_2(:, 6) = matmul(m12+transpose(m12), q_2(1:24))

                bbartrans_a(:, 1) = matmul(m11               , q_a(1:24))
                bbartrans_a(:, 2) = matmul(m22               , q_a(1:24))
                bbartrans_a(:, 3) = matmul(m33               , q_a(1:24))
                bbartrans_a(:, 4) = matmul(m23+transpose(m23), q_a(1:24))
                bbartrans_a(:, 5) = matmul(m13+transpose(m13), q_a(1:24))
                bbartrans_a(:, 6) = matmul(m12+transpose(m12), q_a(1:24))
             
                bbar_2 = transpose(bbartrans_2)

                bbar_a = transpose(bbartrans_a)

                bhat(1, 1) = xi(1)             
                bhat(2, 2) = xi(2)
                bhat(3, 5) = xi(3)             
                bhat(3, 6) = xi(3)*xi(1)
                bhat(3, 7) = xi(3)*xi(2)
                bhat(3, 8) = xi(3)*xi(1)*xi(2)
                bhat(6, 3) = xi(1)
                bhat(6, 4) = xi(2)
             
                bhat = this%sqrtgm_0/sqrtg_0*matmul(curva2curvb6matrix(curva2curvb3matrix(this%gmsub1_0, this%gmsub2_0, this%gmsub3_0, gsup1_0, gsup2_0, gsup3_0)), bhat)
             
                btotal_a(:,  1:24) = bbar_a
                btotal_a(:, 25:32) = bhat

                btotal_2(:,  1:24) = bbar_2
                btotal_2(:, 25:32) = bhat

! TODO strain calculations might be optimized
                strainc_1(1) = 0.5d0*(dot_product(q_1(1:24), matmul(m11, q_1(1:24)))-dot_product(this%q_0(1:24), matmul(m11, this%q_0(1:24))))
                strainc_1(2) = 0.5d0*(dot_product(q_1(1:24), matmul(m22, q_1(1:24)))-dot_product(this%q_0(1:24), matmul(m22, this%q_0(1:24))))
                strainc_1(3) = 0.5d0*(dot_product(q_1(1:24), matmul(m33, q_1(1:24)))-dot_product(this%q_0(1:24), matmul(m33, this%q_0(1:24))))
                strainc_1(4) =       (dot_product(q_1(1:24), matmul(m23, q_1(1:24)))-dot_product(this%q_0(1:24), matmul(m23, this%q_0(1:24))))
                strainc_1(5) =       (dot_product(q_1(1:24), matmul(m13, q_1(1:24)))-dot_product(this%q_0(1:24), matmul(m13, this%q_0(1:24))))
                strainc_1(6) =       (dot_product(q_1(1:24), matmul(m12, q_1(1:24)))-dot_product(this%q_0(1:24), matmul(m12, this%q_0(1:24))))
                
                strainc_2(1) = 0.5d0*(dot_product(q_2(1:24), matmul(m11, q_2(1:24)))-dot_product(this%q_0(1:24), matmul(m11, this%q_0(1:24))))
                strainc_2(2) = 0.5d0*(dot_product(q_2(1:24), matmul(m22, q_2(1:24)))-dot_product(this%q_0(1:24), matmul(m22, this%q_0(1:24))))
                strainc_2(3) = 0.5d0*(dot_product(q_2(1:24), matmul(m33, q_2(1:24)))-dot_product(this%q_0(1:24), matmul(m33, this%q_0(1:24))))
                strainc_2(4) =       (dot_product(q_2(1:24), matmul(m23, q_2(1:24)))-dot_product(this%q_0(1:24), matmul(m23, this%q_0(1:24))))
                strainc_2(5) =       (dot_product(q_2(1:24), matmul(m13, q_2(1:24)))-dot_product(this%q_0(1:24), matmul(m13, this%q_0(1:24))))
                strainc_2(6) =       (dot_product(q_2(1:24), matmul(m12, q_2(1:24)))-dot_product(this%q_0(1:24), matmul(m12, this%q_0(1:24))))
                
                straine_1 = matmul(bhat, q_1(25:32))
                
                straine_2 = matmul(bhat, q_2(25:32))
                
                strain_1 = strainc_1+straine_1 ! addition of compatible and enhanced strains.
                
                strain_2 = strainc_2+straine_2 ! addition of compatible and enhanced strains.
                
                curv2cart3 = curva2curvb3matrix(this%e1_0, this%e2_0, this%e3_0, gsup1_0, gsup2_0, gsup3_0)
                
                curv2cart6 = curva2curvb6matrix(curv2cart3)
                
                cscriptcurv = matmul(transpose(curv2cart6), matmul(this%laminate%plys(l)%cscript, curv2cart6))

                deltastrain = strain_2-strain_1

                ! Initialisierung der Dissipationsanteile
                this%dissipation%sdiss      = 0.0d0
                this%dissipation%vdiss      = 0.0d0
                this%dissipation%cscript_vv = 0.0d0
                this%dissipation%cscript_ev = 0.0d0
                this%dissipation%cscript_ve = 0.0d0
                this%dissipation%cscript_ee = 0.0d0                
                
                if (this%simutype == 0 .or. this%simutype == 3) then
                  call this%dissipation%evaluate(strain_1, strain_2, vh_1, vh_2, cscriptcurv, this%laminate%plys(l)%density, this%laminate%plys(l)%alpha_s, this%laminate%plys(l)%alpha_v, this%deltat)
                end if 
             
                stress_1 = matmul(cscriptcurv, strain_1)
                
                stress_2 = matmul(cscriptcurv, strain_2)
                
                stress_a = 0.5d0*(stress_1+stress_2)+this%dissipation%sdiss
                
                btotal_a_t = transpose(btotal_a)

! TODO calculation of w1 and w2 take lot of time

                w1 = matmul(btotal_a_t, matmul(cscriptcurv+2.0d0*this%dissipation%cscript_ee, btotal_2))
                
                w2(1:24, 1:24) = stress_a(1)* m11+&
                                 stress_a(2)* m22+&
                                 stress_a(3)* m33+&
                                 stress_a(4)*(m23+transpose(m23))+&
                                 stress_a(5)*(m13+transpose(m13))+&
                                 stress_a(6)*(m12+transpose(m12))
             
                this%penergy = this%penergy+&
                               wijk_times_sqrtg_0*0.5d0*dot_product(strain_2, stress_2)

                this%denergy = this%denergy+&
                               wijk_times_sqrtg_0*this%dissipation%energy
             
                nqt = transpose(nq)
                
                this%fqdyn   = this%fqdyn     +&
                               wijk_times_sqrtg_0*(this%laminate%plys(l)%density*matmul(nqt, (vh_2-vh_1)/this%deltat))

                this%fqint   = this%fqint     +&
                               wijk_times_sqrtg_0*(matmul(btotal_a_t, stress_a))
                
                this%fv      = this%fv     +&
                               wijk_times_sqrtg_0*(this%laminate%plys(l)%density*matmul(nqt, (qh_2-qh_1)/this%deltat-(vh_a+this%dissipation%vdiss)))
                
! geometrical and material stiffness matrices for linear buckling analysis
                kqq_mat = 0.0d0
                kqq_geo = 0.0d0
                
                if (this%flag_kgeo_on .eqv. .TRUE.) then
                  kqq_geo = wijk_times_sqrtg_0*w2*0.5d0
                end if
                
                if (this%flag_kmat_on .eqv. .TRUE.) then
                  kqq_mat = wijk_times_sqrtg_0*w1*0.5d0
                end if
                
                this%kqq_geo = this%kqq_geo + kqq_geo
                this%kqq_mat = this%kqq_mat + kqq_mat
                nq_square_dt = matmul(nqt, nq)/this%deltat
                this%kqv     = this%kqv    +&
                               wijk_times_sqrtg_0*(this%laminate%plys(l)%density* nq_square_dt                     +matmul(btotal_a_t         , matmul(this%dissipation%cscript_ev, nq      )))
                this%kvv     = this%kvv    +&
                               wijk_times_sqrtg_0* this%laminate%plys(l)%density*(nq_square_dt*(-0.5d0)* this%deltat -matmul(nqt                , matmul(this%dissipation%cscript_vv, nq      )))
                this%kvq     = this%kvq    +&
                               wijk_times_sqrtg_0* this%laminate%plys(l)%density*(nq_square_dt                     -matmul(nqt                , matmul(this%dissipation%cscript_ve, btotal_2)))

             end do
             
          end do
          
       end do
       
    end do
    
    this%fv (25:32)        = v_2(25:32)
    this%kvv(25:32, 25:32) = eye(8)
   
! For modal and buckling analysis: static condensation of strain terms: kqq* = kqq - B*A_eq = kqq - kqe * kee^-1 * keq = kqq - kqe*(A_eq); kee*A_eq = keq
    if ( (this%simutype == 2) .or. (this%simutype == 3) ) then                                          ! modal, buckling
      this%kqq = this%kqq_geo + this%kqq_mat
      kqq_temp = 0.0d0
      kqq_temp = this%kqq
      kee(:,:) = kqq_temp(25:32,25:32)                                    
      RHS(:,:) = kqq_temp(25:32,1:24)                                     
      allocate(ipiv(8))
      call dgesv(8, 24, kee, 8, ipiv, RHS, 8, info)                       ! Solving system of equation
      kqq_sc = kqq_temp(1:24,1:24) - matmul(kqq_temp(1:24,25:32),RHS)     ! kqq - kqe*A_eq
      this%kqq(1:24,1:24) = kqq_sc                                        ! overwriting new kqq including condensation at position in kqq-components
      deallocate(ipiv)
    else                                                                  ! Static and dynamic analysis
      this%kqq = this%kqq_geo + this%kqq_mat
    end if
    
return
    
  end subroutine shellinternalterms

  function nphimatrix(n) result(nphi)

    implicit none

    real(kind = 8), intent(in) :: n(4)

    real(kind = 8) :: nphi(3, 24)

    nphi(:, :) = 0.0d0

    nphi(:,  1: 3) = i_33*n(1) 
    nphi(:,  7: 9) = i_33*n(2)
    nphi(:, 13:15) = i_33*n(3)
    nphi(:, 19:21) = i_33*n(4)
           
    return
    
  end function nphimatrix

  function ndirmatrix(n) result(ndir)

    implicit none

    real(kind = 8), intent(in) :: n(4)

    real(kind = 8) :: ndir(3, 24)

    ndir(:, :) = 0.0d0

    ndir(:,  4: 6) = i_33*n(1)
    ndir(:, 10:12) = i_33*n(2)
    ndir(:, 16:18) = i_33*n(3)
    ndir(:, 22:24) = i_33*n(4)

    return
    
  end function ndirmatrix

  subroutine unitsquareshapefunctions(n, d1n, d2n, xi)

    implicit none

    real(kind = 8), intent(in) :: xi(:)

    real(kind = 8), dimension(4), intent(out) :: n, d1n, d2n
    
    n(1) = 0.25d0*(1+xi(1))*(1+xi(2))
    n(2) = 0.25d0*(1-xi(1))*(1+xi(2))
    n(3) = 0.25d0*(1-xi(1))*(1-xi(2))
    n(4) = 0.25d0*(1+xi(1))*(1-xi(2))
    
    d1n(1) = 0.25d0*(1+xi(2))
    d1n(2) =-0.25d0*(1+xi(2))
    d1n(3) =-0.25d0*(1-xi(2))
    d1n(4) = 0.25d0*(1-xi(2))   

    d2n(1) = 0.25d0*(1+xi(1))
    d2n(2) = 0.25d0*(1-xi(1))
    d2n(3) =-0.25d0*(1-xi(1))
    d2n(4) =-0.25d0*(1+xi(1))
    
    return
    
  end subroutine unitsquareshapefunctions

  subroutine mmatrices(m11, m22, m33, m23, m13, m12, xi, zeta, d3zeta, nphi, ndir, d1nphi, d1ndir, d2nphi, d2ndir)

    implicit none

    real(kind = 8), intent(in) :: xi(:), zeta, d3zeta

    real(kind = 8), dimension(3, 24), intent(in) :: nphi, ndir, d1nphi, d1ndir, d2nphi, d2ndir
    real(kind = 8), dimension(24, 3) :: d1nphi_t, d1ndir_t, d2nphi_t, d2ndir_t

    real(kind = 8), dimension(24, 24), intent(out) :: m11, m22, m33, m23, m13, m12

    real(kind = 8), dimension(3, 24) :: ndira, ndirb
    real(kind = 8), dimension(3, 24) :: ndirc, ndird
    
    real(kind = 8), dimension(3, 24) :: d1nphia, d1nphic
    real(kind = 8), dimension(3, 24) :: d2nphib, d2nphid
        
    real(kind = 8), dimension(24, 24) :: m23zth, m23fst, m13zth, m13fst

! TODO these terms are invariant. initialize them in the module
    ndira(:, :) = 0.0d0
    ndira(1:3,  4: 6) = 0.5d0*i_33
    ndira(1:3, 10:12) = 0.5d0*i_33

    ndirb(:, :) = 0.0d0
    ndirb(1:3, 10:12) = 0.5d0*i_33
    ndirb(1:3, 16:18) = 0.5d0*i_33

    ndirc(:, :) = 0.0d0
    ndirc(1:3, 16:18) = 0.5d0*i_33
    ndirc(1:3, 22:24) = 0.5d0*i_33

    ndird(:, :) = 0.0d0
    ndird(1:3,  4: 6) = 0.5d0*i_33
    ndird(1:3, 22:24) = 0.5d0*i_33

    d1nphia(:, :) = 0.0d0
    d1nphia(1:3,  1: 3) = 0.5d0*i_33
    d1nphia(1:3,  7: 9) =-0.5d0*i_33

    d2nphib(:, :) = 0.0d0
    d2nphib(1:3,  7: 9) = 0.5d0*i_33
    d2nphib(1:3, 13:15) =-0.5d0*i_33

    d1nphic(:, :) = 0.0d0
    d1nphic(1:3, 13:15) =-0.5d0*i_33
    d1nphic(1:3, 19:21) = 0.5d0*i_33

    d2nphid(:, :) = 0.0d0
    d2nphid(1:3,  1: 3) = 0.5d0*i_33
    d2nphid(1:3, 19:21) =-0.5d0*i_33
    
    d1nphi_t = transpose(d1nphi)
    d1ndir_t = transpose(d1ndir)
    d2nphi_t = transpose(d2nphi)
    d2ndir_t = transpose(d2ndir)
    
    m11 =          matmul(d1nphi_t, d1nphi)  &                                    
         +zeta   *(matmul(d1nphi_t, d1ndir)+ &
                   matmul(d1ndir_t, d1nphi)) & 
         +zeta**2* matmul(d1ndir_t, d1ndir)

    m22 =          matmul(d2nphi_t, d2nphi)  &
         +zeta   *(matmul(d2nphi_t, d2ndir)+ &
                   matmul(d2ndir_t, d2nphi)) &
         +zeta**2* matmul(d2ndir_t, d2ndir)

    m33(:, :) = 0.0d0

    m33( 4: 6,  4: 6) = i_33*ndir(1,  4)
    m33(10:12, 10:12) = i_33*ndir(1, 10)
    m33(16:18, 16:18) = i_33*ndir(1, 16)
    m33(22:24, 22:24) = i_33*ndir(1, 22)
    
    m33 = m33*d3zeta**2
    
    m23zth = 0.5d0*((1-xi(1))*matmul(transpose(d2nphib), ndirb)+&
                    (1+xi(1))*matmul(transpose(d2nphid), ndird))! ZeroTH-order
    
    m23fst = matmul(d2ndir_t, ndir)! FirST-order

    m23 = (m23zth+zeta*m23fst)*d3zeta

    m13zth = 0.5d0*((1-xi(2))*matmul(transpose(d1nphic), ndirc)+&
                    (1+xi(2))*matmul(transpose(d1nphia), ndira))! ZeroTH-order
        
    m13fst = matmul(d1ndir_t, ndir)! FirST-order

    m13 = (m13zth+zeta*m13fst)*d3zeta

    m12 =          matmul(d1nphi_t, d2nphi)  &                                    
         +zeta   *(matmul(d1nphi_t, d2ndir)+ &
                   matmul(d1ndir_t, d2nphi)) & 
         +zeta**2* matmul(d1ndir_t, d2ndir)
    
    return

  end subroutine mmatrices
        
  function shellmassmatrix(this) result(massmatrix)

    implicit none
    
    class(shell_element), intent(inout) :: this
    
    real(kind = 8) :: massmatrix(24, 24)

    real(kind = 8) :: sqrtg_0, volume_0
    real(kind = 8), dimension(3) :: gsub1_0, gsub2_0, gsub3_0

    real(kind = 8), dimension(3, 24) :: nphi, ndir, d1nphi, d1ndir, d2nphi, d2ndir

    real(kind = 8) :: zeta, d3zeta
    
    integer :: ngaussianpoints = 2
    real(kind = 8), dimension(2) :: xigauss, wigauss
    real(kind = 8) :: xi(3)
    real(kind = 8), dimension(4) :: n, d1n, d2n

    real(kind = 8) :: nq(3, 24)
    
    integer :: i, j, k, l

    real(kind = 8) :: wijk_times_sqrtg_0
                
    xigauss = 0.577350269189626d0*[-1.0d0, 1.0d0]

    wigauss =                     [ 1.0d0, 1.0d0]
    
    massmatrix(:, :) = 0.0d0
    
    volume_0 = 0.0d0

    do l = 1, this%laminate%nplys
     
        d3zeta = this%laminate%thickness*0.5d0
        
       do k = 1, ngaussianpoints

          do j = 1, ngaussianpoints

             do i = 1, ngaussianpoints

                xi = [xigauss(i), xigauss(j), xigauss(k)]

                zeta = this%laminate%boundaries(l)*(1.0d0-xi(3))*0.5d0+this%laminate%boundaries(l+1)*(1.0d0+xi(3))*0.5d0
                
                call unitsquareshapefunctions(n, d1n, d2n, xi)
                
                nphi   = nphimatrix(n)
                d1nphi = nphimatrix(d1n)
                d2nphi = nphimatrix(d2n)
                
                ndir   = ndirmatrix(n)
                d1ndir = ndirmatrix(d1n)
                d2ndir = ndirmatrix(d2n)             
                
                nq = nphi+zeta*ndir
                
                gsub1_0 = matmul(d1nphi+zeta*d1ndir, this%q_0(1:24))
                
                gsub2_0 = matmul(d2nphi+zeta*d2ndir, this%q_0(1:24))
                
                gsub3_0 = d3zeta*matmul(ndir, this%q_0(1:24))
                
                sqrtg_0 = dot_product(cross(gsub1_0, gsub2_0), gsub3_0)
                                
                wijk_times_sqrtg_0 = wigauss(i)*wigauss(j)*wigauss(k)*sqrtg_0*this%laminate%plys(l)%thickness/this%laminate%thickness
                
                massmatrix = massmatrix+wijk_times_sqrtg_0*this%laminate%plys(l)%density*matmul(transpose(nq), nq)
                
                volume_0   = volume_0  +wijk_times_sqrtg_0             
                
             end do
             
          end do
          
       end do

    end do

    this%length_0 = dsqrt(volume_0/this%laminate%thickness)! volume_0**(1.0d0/3.0d0)
    
    this%mass_matrix = massmatrix
    return
    
  end function shellmassmatrix
  
end module class_shell_element
