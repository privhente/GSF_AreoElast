module class_inflow_aero
  use my_FileIO_functions
  implicit none

  type :: air_flow
      integer :: ifileType, NumGrid_Y, NumGrid_Z, NumOutSteps
      real(kind = 8) :: gridheight, gridwidth
      real(kind = 8) :: dir(3), vinfinity(3), center(3)
      real(kind = 8) :: density, intensity, duration, totalt, deltat, vmean
      real(kind = 8), allocatable :: arrvinf(:,:,:,:)
      character(len = :), allocatable :: sort
      character(len = :), allocatable :: strfilename_windfield
      logical :: boolean_external  = .FALSE.
      logical :: boolean_actualize = .TRUE.
      logical :: boolean_abort = .FALSE.
    contains
      procedure :: airflowvinfinity
      procedure :: getairflowDatafromfile
      procedure :: ReadFullFieldBinary_Bladed_wnd
      procedure :: writeWindFields_ASCII
    end type air_flow
  
  contains
  
!< subroutine to read air flow data from file, depending on the data format (.wnd, .bts, ...)
  subroutine getairflowDatafromfile(this)
    implicit none
    
    class(air_flow) :: this
    character(len = :), allocatable :: str_ext
    integer :: ifileunit, io_error, inz_ext
    
    allocate(str_ext, source = '-')
    
    !< finding out the character of file
    inz_ext = index(this%strfilename_windfield, '.wnd')
    if (inz_ext .ne. 0) then 
      deallocate(str_ext)
      allocate(str_ext, source = this%strfilename_windfield(inz_ext:len(this%strfilename_windfield))) 
    end if
    
    inz_ext = index(this%strfilename_windfield, '.bts')
    if (inz_ext .ne. 0) then 
      deallocate(str_ext)
      allocate(str_ext, source = this%strfilename_windfield(inz_ext:len(this%strfilename_windfield))) 
    end if
  
    select case (str_ext)
      case('.wnd')
        !< open binary file
        call GetNewUnit(ifileunit)
        call OpenBOutFileToRead(ifileunit, this%strfilename_windfield, io_error)
        if (io_error .eq. 0) then
          print*, 'Parsing air flow file, opening ', this%strfilename_windfield
          call this%ReadFullFieldBinary_Bladed_wnd(ifileunit)
          if (this%boolean_abort) then
            print*, 'Error: while in reading air flow file ', this%strfilename_windfield
            return
          else
            print '(A17,1X,I5)', 'grid size y:', this%NumGrid_Y
            print '(A17,1X,I5)', 'grid size z:', this%NumGrid_Z
            print '(A17,1X,F10.5)', 'grid height:', this%gridheight
            print '(A17,1X,F10.5)', 'grid width:', this%gridwidth
            print '(A17,1X,F10.5)', 'data time:', this%totalt
            print '(A17,1X,F10.5)', 'delta t:', this%deltat
            print '(A17,1X,F10.5)', 'vmean:', this%vmean
          end if
          this%boolean_external = .TRUE.
        else
          this%boolean_abort = .TRUE.
          print*, 'Error in reading air flow file ', this%strfilename_windfield
        end if
        close(ifileunit)

      case('.bts')
        !< open binary file
        call GetNewUnit(ifileunit)
        !call OpenBOutFileToRead(ifileunit, this%strfilename_windfield, io_error)
        io_error = -1
        if (io_error .eq. 0) then
          print*, 'Parsing air flow file, opening ', this%strfilename_windfield
        else
          this%boolean_abort = .TRUE.
          print*, 'Error in reading air flow file ', this%strfilename_windfield
        end if
        close(ifileunit)
      case default
        this%boolean_abort = .TRUE.
        print*, 'Error in reading air flow file ', this%strfilename_windfield
        return
      end select 
    
    ! write wind fields to file in ascii-format
    !if (allocated(this%arrVinf)) call this%writeWindFields_ASCII()
    
    return
  end subroutine getairflowDatafromfile

  !< subroutine to write read wind fields from binary to read-able ascii-format
  subroutine writeWindFields_ASCII(this)
    implicit none
    
    class(air_flow) :: this
    integer :: it, iz, iy, ifileunit_u1, ifileunit_u2, ifileunit_u3
    character(len=7) :: str_status = 'replace'
    character(len=6) :: str_position = 'rewind'

1000 format(1000000f10.3)
     
    call GetNewUnit(ifileunit_u1)
    open(unit = ifileunit_u1, file = this%strfilename_windfield // '_u1.dres', action = 'write', status = str_status, POSITION = str_position)
    
    call GetNewUnit(ifileunit_u2)
    open(unit = ifileunit_u2, file = this%strfilename_windfield // '_u2.dres', action = 'write', status = str_status, POSITION = str_position)
    
    call GetNewUnit(ifileunit_u3)
    open(unit = ifileunit_u3, file = this%strfilename_windfield // '_u3.dres', action = 'write', status = str_status, POSITION = str_position)
    
    write(ifileunit_u1, '(A35)') '!! Summary of wind data information'
    write(ifileunit_u1, '(A17,1X,I5)') 'grid size y:', this%NumGrid_Y
    write(ifileunit_u1, '(A17,1X,I5)') 'grid size z:', this%NumGrid_Z
    write(ifileunit_u1, '(A17,2X,F10.5)') 'grid height:', this%gridheight
    write(ifileunit_u1, '(A16,3X,F10.5)') 'grid width:', this%gridwidth
    write(ifileunit_u1, '(A15,4X,F10.5)') 'data time:', this%totalt
    write(ifileunit_u1, '(A13,6X,F10.5)') 'delta t:', this%deltat
    write(ifileunit_u1, '(A11,8X,F10.5)') 'vmean:', this%vmean
    
    write(ifileunit_u2, '(A35)') '!! Summary of wind data information'
    write(ifileunit_u2, '(A17,1X,I5)') 'grid size y:', this%NumGrid_Y
    write(ifileunit_u2, '(A17,1X,I5)') 'grid size z:', this%NumGrid_Z
    write(ifileunit_u2, '(A17,2X,F10.5)') 'grid height:', this%gridheight
    write(ifileunit_u2, '(A16,3X,F10.5)') 'grid width:', this%gridwidth
    write(ifileunit_u2, '(A15,4X,F10.5)') 'data time:', this%totalt
    write(ifileunit_u2, '(A13,6X,F10.5)') 'delta t:', this%deltat
    write(ifileunit_u2, '(A11,8X,F10.5)') 'vmean:', this%vmean
    
    write(ifileunit_u3, '(A35)') '!! Summary of wind data information'
    write(ifileunit_u3, '(A17,1X,I5)') 'grid size y:', this%NumGrid_Y
    write(ifileunit_u3, '(A17,1X,I5)') 'grid size z:', this%NumGrid_Z
    write(ifileunit_u3, '(A17,2X,F10.5)') 'grid height:', this%gridheight
    write(ifileunit_u3, '(A16,3X,F10.5)') 'grid width:', this%gridwidth
    write(ifileunit_u3, '(A15,4X,F10.5)') 'data time:', this%totalt
    write(ifileunit_u3, '(A13,6X,F10.5)') 'delta t:', this%deltat
    write(ifileunit_u3, '(A11,8X,F10.5)') 'vmean:', this%vmean

    !< loop over time steps
    do it = 1, this%NumOutSteps                     ! Use only the number of timesteps requested originally
      write(ifileunit_u1, '(A8,1X,F10.5)') '!! time:', this%deltat*(it-1)
      write(ifileunit_u2, '(A8,1X,F10.5)') '!! time:', this%deltat*(it-1)
      write(ifileunit_u3, '(A8,1X,F10.5)') '!! time:', this%deltat*(it-1)
      do iz = 1, this%NumGrid_Z
        do iy = 1, this%NumGrid_Y
          write(ifileunit_u1, 1000, advance='no') this%arrVinf(iy,iz,it,1) ! u1 - in x-dierction
          write(ifileunit_u2, 1000, advance='no') this%arrVinf(iy,iz,it,2) ! u2 - in y-dierction
          write(ifileunit_u3, 1000, advance='no') this%arrVinf(iy,iz,it,3) ! u3 - in z-dierction
        end do
        write(ifileunit_u1, 1000, advance='yes')
        write(ifileunit_u2, 1000, advance='yes')
        write(ifileunit_u3, 1000, advance='yes')
      end do
    end do
    
    close(ifileunit_u1)
    close(ifileunit_u2)
    close(ifileunit_u3)
  
    return
  end subroutine
  
  
  !< read full-field binary file in Bladed (.wnd) file format
  subroutine ReadFullFieldBinary_Bladed_wnd(this, ifileunit)
    implicit none
    
    class(air_flow) :: this
    integer, intent(in) :: ifileunit
    
    !< parameter to read from file
    integer(kind = B2Ki) :: blade_format, improved_von_karman
    integer(kind = B4Ki) :: number_wind_components
    real(kind = SiKi) :: latitude, roughness_length, reference_height, TI(3)
    real(kind = SiKi) :: GridRes_Z, GridRes_Y, UHub
    
    integer(kind = B4Ki) :: NumOutStepshalf
    real( kind = SiKi) :: mean_Vhh, temp1, temp2, temp3
    integer(kind = B4Ki) :: temp4, RandSeed
    integer(kind = B4Ki) :: temp5, temp6, temp7, temp8, temp9, temp10
    integer(kind = B2Ki), allocatable :: TmpVarray(:)
    real(kind = SiKi) :: U_C1, U_C2, V_C, W_C
    
    !real(kind = 8), allocatable :: arrGrid_Y(:), arrGrid_Z(:)
    integer :: it, ip, iy, iz, r1, indicesv(3)
    integer :: io_error
    integer(kind = B4Ki) :: NumGrid_Y, NumGrid_Z
    
    !< reading the header
    read(ifileunit, iostat = io_error) blade_format                    ! -99 = New Bladed format
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) improved_von_karman             ! 4 = improved von karman (but needed for next 7 inputs)
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) number_wind_components          ! 3 = number of wind components
    if (io_error .ne. 0) goto 10

    read(ifileunit, iostat = io_error) latitude                        ! Latitude (degrees)
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) roughness_length                ! Roughness length (m)
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) reference_height                ! Reference Height (m) ( Z(1) + GridHeight / 2.0 ) !This is the vertical center of the grid
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) TI(1)                           ! Longitudinal turbulence intensity (%)
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) TI(2)                           ! Lateral turbulence intensity (%)
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) TI(3)                           ! Vertical turbulence intensity (%)
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) GridRes_Z                       ! grid spacing in vertical direction, in m
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) GridRes_Y                       ! grid spacing in lateral direction, in m
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) UHub                            ! grid spacing in longitudinal direciton, in m
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) NumOutStepshalf                 ! half the number of points in along wind direction
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) mean_Vhh                        ! the mean wind speed in m/s at hub height
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp1                           ! the vertical length scale of the longitudinal component in m
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp2                           ! the lateral length scale of the longitudinal component in m
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp3                           ! the longitudinal length scale of the longitudinal component in m
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp4                           ! an unused integer
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) RandSeed                        ! the random number seed
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) NumGrid_Z                       ! the number of grid points vertically
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) NumGrid_Y                       ! the number of grid points laterally
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp5                           ! the vertical length scale of the lateral component, not used
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp6                           ! the lateral length scale of the lateral component, not used
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp7                           ! the longitudinal length scale of the lateral component, not used
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp8                           ! the vertical length scale of the vertical component, not used
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp9                           ! the lateral length scale of the vertical component, not used
    if (io_error .ne. 0) goto 10
    
    read(ifileunit, iostat = io_error) temp10                          ! the longitudinal length scale of the vertical component, not used
    if (io_error .ne. 0) goto 10
    
    !< some numbers for re-normalizing the data
    U_C1 = 1000.0 / (mean_Vhh*TI(1)/100.0)
    U_C2 = 1000.0 / (TI(1)/100)
    V_C  = 1000.0 / (mean_Vhh*TI(2)/100.0)
    W_C  = 1000.0 / (mean_Vhh*TI(3)/100.0)
    
    !< grid data
    this%NumGrid_Y   = NumGrid_Y
    this%NumGrid_Z   = NumGrid_Z
    this%gridheight  = GridRes_Z * DBLE((this%NumGrid_Z - 1))
    this%gridwidth   = GridRes_Y * DBLE((this%NumGrid_Y - 1))
    this%deltat      = UHub/mean_Vhh
    this%vmean       = mean_Vhh
    this%NumOutSteps = 2*NumOutStepshalf
    this%totalt      = this%deltat*DBLE(this%NumOutSteps)
    
    ! allocating velocitiy array
    allocate(TmpVarray(3*this%NumGrid_Y*this%NumGrid_Z))
    allocate(this%arrVinf(this%NumGrid_Y,this%NumGrid_Z,this%NumOutSteps,3))
    
    !< loop over time steps
    do it = 1, this%NumOutSteps                     ! Use only the number of timesteps requested originally
      read(ifileunit) TmpVarray                     ! the longitudinal length scale of the vertical component, not used
      ip = 1
      do iz = 1, this%NumGrid_Z
        do iy = 1, this%NumGrid_Y
          indicesv = (/ (r1, r1 = 3*(ip-1)+1, 3*(ip-1)+3) /)
          !< re-calculation of normalized velocity components according to "TS_FileIO.f90" routine "WrBinBLADED" in TurbSim source code
          this%arrVinf(iy,iz,it,1) = (real(TmpVarray(indicesv(1))) + U_C2 ) / U_C1 ! u - in x-dierction
          this%arrVinf(iy,iz,it,2) =  real(TmpVarray(indicesv(2))) / V_C           ! v - in y-dierction
          this%arrVinf(iy,iz,it,3) =  real(TmpVarray(indicesv(3))) / W_C           ! w - in z-dierction
          ip = ip + 1
        end do
      end do
    end do
    
    return
10  continue
    this%boolean_abort = .TRUE.
    return
  end subroutine ReadFullFieldBinary_Bladed_wnd
  
!< function to get external flow field for a point with coordinates q
  function airflowvinfinity(this, q, time) result(airflow_vinfinity)
    implicit none
    class(air_flow), intent(in) :: this
    real(kind = 8), intent(in) :: q(3)
    real(kind = 8), intent(in) :: time
    integer :: igridID_y, igridID_z, igridID_t
    real(kind = 8) :: temp_time
    real(kind = 8) :: airflow_vinfinity(3) 
    real(kind = 8) :: gridID_y, gridID_z, gridID_t
    real(kind = 8), dimension(3) :: vy1_z1_t1, vy2_z1_t1, vy1_z2_t1, vy2_z2_t1, vy1_z1_t2, vy2_z1_t2, vy1_z2_t2, vy2_z2_t2, v_t1, v_t2
    real(kind = 8) :: fy, fz, ft

    temp_time = time
    if (time .gt. this%totalt) temp_time = this%totalt
  
    !< take velocity components from array according to indices
    !< if point is outside of grid then assign constant velocity according to simulation settings, else take from external free-field velocity field
    if (abs(q(2)) .gt. (this%gridwidth*0.5d0 + this%center(2)) .or. abs(q(3)) .gt. (this%gridheight*0.5d0 + this%center(3))) then
      airflow_vinfinity = this%vinfinity
    else
      !< determine indices for grid and time depend velocity array for point and time step
      gridID_y = DBLE(this%NumGrid_Y - 1) / this%gridwidth  * (q(2) - this%center(2) + this%gridwidth*0.5d0) + 1.0d0
      gridID_z = DBLE(this%NumGrid_Z - 1) / this%gridheight * (q(3) - this%center(3) + this%gridheight*0.5d0) + 1.0d0
      gridID_t = DBLE(this%NumOutSteps - 1) / this%totalt * temp_time + 1.0d0
      
      igridID_y  = int(gridID_y)
      igridID_z  = int(gridID_z)
      igridID_t  = int(gridID_t)
      
      !< factor for interpolation
      fy = abs(gridID_y - dble(igridID_y))
      fz = abs(gridID_z - dble(igridID_z))
      ft = abs(gridID_t - dble(igridID_t))
      
      !< factor for interpolation
      vy1_z1_t1 = this%arrVinf(igridID_y  ,igridID_z  ,igridID_t,1:3)
      vy2_z1_t1 = this%arrVinf(igridID_y+1,igridID_z  ,igridID_t,1:3)
      vy1_z2_t1 = this%arrVinf(igridID_y  ,igridID_z+1,igridID_t,1:3)
      vy2_z2_t1 = this%arrVinf(igridID_y+1,igridID_z+1,igridID_t,1:3)
      
      !< interpolation at the grid for values at time 1
      v_t1(1:3) = fun_interpolationatplane(vy1_z1_t1, vy2_z1_t1, vy1_z2_t1, vy2_z2_t1, fy, fz)
      
      !< interpolation at the grid for values at time 2, if necessary
      v_t2(1:3) = 0.0d0
      if (ft .ne. 0.0d0) then
        vy1_z1_t2 = this%arrVinf(igridID_y  ,igridID_z  ,igridID_t+1,1:3)
        vy2_z1_t2 = this%arrVinf(igridID_y+1,igridID_z  ,igridID_t+1,1:3)
        vy1_z2_t2 = this%arrVinf(igridID_y  ,igridID_z+1,igridID_t+1,1:3)
        vy2_z2_t2 = this%arrVinf(igridID_y+1,igridID_z+1,igridID_t+1,1:3)
        v_t2(1:3) = fun_interpolationatplane(vy1_z1_t2, vy2_z1_t2, vy1_z2_t2, vy2_z2_t2, fy, fz)
      end if

      !< resultant airflow velocity at point
      airflow_vinfinity(1:3) = (v_t2 - v_t1)*ft + v_t1
    end if
    
    return
  end function airflowvinfinity
  
  !< function for interpolating between four points considering
  !< vp(y) = (v2-v1)*ratio+v1, with ratio = (yp-y1)/(y2-y1)
  function fun_interpolationatplane(vy1_z1, vy2_z1, vy1_z2, vy2_z2, fy, fz) result(vp)
    implicit none
    real(kind = 8), intent(in), dimension(3) ::  vy1_z1, vy2_z1, vy1_z2, vy2_z2
    real(kind = 8) :: fy, fz
    real(kind = 8) :: vp(3)
    
    vp = (vy1_z1 + (vy2_z1 - vy1_z1)*fy) + ( (vy1_z2 + (vy2_z2 - vy1_z2)*fy) - (vy1_z1 + (vy2_z1 - vy1_z1)*fy) )*fz
    
    return
  end function fun_interpolationatplane
  
end module class_inflow_aero
