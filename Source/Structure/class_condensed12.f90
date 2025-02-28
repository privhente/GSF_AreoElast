module class_condensed12
  use my_constants_structure, only: i_33, e1, e2, e3

  implicit none

  type :: condensed_12_curve                                                   !< sub user defined variable for condensed element type: property
    real(kind = 8), allocatable :: data_points(:,:)                            !< array of load deflection diagram
    integer :: npoints                                                         !< number of points in load deflection diagram
  end type condensed_12_curve

  type :: condensed_12
    integer :: node                                                             !< node indices for condensed element
    integer :: propertynumber                                                   !< load deflection diagram
    real(kind = 8) :: deltat
    real(kind = 8) :: dir(3)                                                    !< direction vector
    real(kind = 8) :: penergy                                                   !< potential energy
    real(kind = 8) :: load, stiffness, xf                                       !< result variables load and stiffness due to load deflection diagram
    integer :: npoints
    real(kind = 8), dimension(12,12) :: kqq, kqv, kvq, kvv
    real(kind = 8), dimension(12) :: fqint, fv
    character(:), allocatable :: sort                                               !< kind of selected condensed element
    real(kind = 8), allocatable :: data_points(:,:), ranges(:,:), const_stiff(:)    !< settings and variables for function linearstiffness
    contains
      procedure :: linearstiffness
  end type condensed_12

contains

! Function to determine linear spring12 stiffness and interpolated spring load
  subroutine linearstiffness(this,xi)
    implicit none

    class(condensed_12), intent(inout) :: this
    real(kind = 8), intent(in) :: xi

    real(kind = 8) :: max_points, min_points
    real(kind = 8) :: z, x2, x1, Fr, Fl, Fm
    integer :: k, np, nranges

    np      = size(this%data_points,1)
    nranges = size(this%ranges,1)

    max_points  = maxval(this%data_points(:,1),1)
    min_points  = minval(this%data_points(:,1),1)

    this%stiffness = 0.0d0
    this%load      = 0.0d0
    this%xf        = xi

    !if (xi .eq. 0.0d0) then
    !  return
    !end if

    if (isnan(xi) .eqv. .false.) then
      if (xi >= max_points) then
        this%stiffness = this%const_stiff(np)
        this%load      = this%stiffness*(xi-this%data_points(np,1))+this%data_points(np,2)
      elseif (xi < min_points) then
        this%stiffness = 0.0d0
        this%load      = 0.0d0
      else
        do k = 1,nranges-1
            x1 = this%ranges(k,1)
            x2 = this%ranges(k+1,1)
            z = (2.0d0*xi-(x1+x2))/(x2-x1)
            if (abs(z) .le. 1.0d0) then
                exit
            end if
        end do
        if (k .eq. 1) then
            Fl = this%data_points(1,2)
            Fr = this%ranges(k+1,2)
            Fm = (Fl+Fr)*0.5d0
            this%stiffness = this%const_stiff(k)
        elseif (k .gt. 1 .and. k .lt. nranges) then
            Fl = this%ranges(k,2)
            Fm = this%data_points(k,2)
            Fr = this%ranges(k+1,2)
            this%stiffness = 0.5d0*(1.0d0-z)*this%const_stiff(k-1) + 0.5d0*(1.0d0+z)*this%const_stiff(k)
        else
            Fl = this%ranges(nranges-1,2)
            Fr = this%ranges(nranges,2)
            Fm = (Fl+Fr)*0.5d0
            this%stiffness = this%const_stiff(nranges-1)
        end if
        this%load = -z*(1.0d0-z)*0.5d0*Fl + (1.0d0-z*z)*Fm + z*(1.0d0+z)*0.5d0*Fr
      end if
    else
      this%load = xi
      this%stiffness = xi ! setting to nan only for initialisation
    end if
    
    return
  end subroutine linearstiffness

end module class_condensed12