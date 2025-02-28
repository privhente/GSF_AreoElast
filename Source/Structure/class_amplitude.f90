module class_amplitude
  use my_constants_structure, only: pi
  implicit none
  
  !< class for force and boundary condition amplitudes
   type :: amplitude
     
     character(len = :), allocatable :: sort !< kind of amplitude selected
     character(len = :), allocatable :: filename !< filename of file to read
     real(kind = 8) :: duration              !< duration of the waveform
     real(kind = 8) :: intensity             !< intensity of the waveform
     
     real(kind = 8), allocatable :: amplitudes(:) !< in case of 'file' this is the amplitude time series
     real(kind = 8), allocatable :: timestamps(:) !< in case of 'file' this is the timestamp series
     
     logical :: boolean_zero = .FALSE.
   contains

   procedure :: amplitude => evaluateamplitude
    
  end type amplitude

contains

  function evaluateamplitude(this, time) result(real_amplitude)
    
    implicit none

    class(amplitude) :: this
    real(kind = 8) :: time
    real(kind = 8) :: real_amplitude
    integer :: io_error, i, nSamples
    
    real_amplitude = 0.0d0
    
    select case (this%sort)
       
    case('pitch')
! dirac delta. Only works for exact times
       if (this%boolean_zero) then
         real_amplitude = 0.0d0
         return
       end if
       if (time .le. this%duration) then
          real_amplitude = 1.0d0
          this%boolean_zero  = .FALSE.
       end if
       
!    case('pitch2')
!! dirac delta. Only works for exact times
!       if (this%boolean_zero) then
!         real_amplitude = 0.0d0
!         return
!       end if
!       if (time .le. this%duration) then
!          real_amplitude = 1.0d0
!          this%boolean_zero  = .FALSE.
!       end if
!       
!    case('pitch3')
!! dirac delta. Only works for exact times
!       if (this%boolean_zero) then
!         real_amplitude = 0.0d0
!         return
!       end if
!       if (time .le. this%duration) then
!          real_amplitude = 1.0d0
!          this%boolean_zero  = .FALSE.
!       end if

    case('delta')
! dirac delta. Only works for exact times
       if (time .eq. this%duration) then
          real_amplitude = 1.0d0
       else
          real_amplitude = 0.0d0
       end if

    case('heaviside')
! step function, normalized
      if( time > this%duration) then
        real_amplitude = 1.0d0
      else
        real_amplitude = 0.0d0
      endif

    case('constant')
! step function, normalized
      if( time >= 0.0d0) then
        real_amplitude = 1.0d0
      else
        real_amplitude = 0.0d0
      endif
      
    case('linear')
! linear increase until duration
       if (time <= this%duration) then
          real_amplitude = time/this%duration	
       else
          real_amplitude = time/this%duration	
       end if
    
    case('linear0')
! linear increase until duration
       if (time <= this%duration) then
          real_amplitude = time/this%duration	
       else
          real_amplitude = 0.0d0
       end if
    
    case('linearC')
! linear increase until duration and then constant
       if (time <= this%duration) then
          real_amplitude = time/this%duration	
       else
          real_amplitude = 1.0d0
       end if
       
    case('triangle')
! linear increase until half of the duration, then linear decay
       if (time < this%duration/2.0d0) then
          real_amplitude = time/(this%duration/2.0d0)	
       else if ((time>=this%duration/2.0d0) .and. (time<this%duration)) then
          real_amplitude = (1.0d0-(time-this%duration/2.0d0)/(this%duration/2.0d0))      
       else
          real_amplitude = 0.0d0
       end if

    case('sinusoidal')
! sine wave
       real_amplitude = dsin(2.0d0*pi*time/this%duration)

    case('cosinusoidal')
! cosine wave
       real_amplitude = dcos(2.0d0*pi*time/this%duration)

    case('square')
! inverted heaviside
       if (time <= this%duration) then
          real_amplitude = 1.0d0
       else if (time >= this%duration) then
          real_amplitude = 0.0d0
       end if
       
    case('random')
! random number
       call random_number(real_amplitude)
       real_amplitude = 2.0d0*real_amplitude - 1.0d0
      
    case('file')
! time series from file
      if(time == 0.0d0) then
        ! initialize time series from file
        print *, "Parsing amplitude file: ",  this%filename
        open(unit = 20000, file = this%filename, status = 'old', iostat = io_error)
        
        if (io_error /= 0) then
          print*, "Error in reading amplitude file. Check filename: ", this%filename, "!"
          stop
        end if
        
        if (io_error == 0) then  ! if the amplitude file exists, then read the information.
            
             ! skip comment lines
             do i = 1, 3
                read(20000, *)
             end do
             
             ! read number of samples
             read(20000, *) nSamples
            
             if (allocated(this%amplitudes)) then
               deallocate(this%amplitudes)
               deallocate(this%timestamps)
             end if
             
             allocate(this%amplitudes(nSamples))
             allocate(this%timestamps(nSamples))
       
             ! skip comment lines
             do i = 1, 3
                read(20000, *)
             end do
       
             ! read time series
             do i = 1, nSamples
                read(20000, *) this%timestamps(i), this%amplitudes(i)
             end do
           
             close(unit = 20000)
           
             ! initialize
             real_amplitude = this%amplitudes(1)
        end if
      else
        ! already initialized
        if (allocated(this%amplitudes)) then
          real_amplitude = this%amplitudes(size(this%amplitudes))
          ! find the time stamp for the current time
          do i = 2, size(this%amplitudes)
            if(this%timestamps(i)>time) then
              ! linear interpolation
              real_amplitude = (this%amplitudes(i) * (time - this%timestamps(i-1)) + this%amplitudes(i-1) * (this%timestamps(i) - time) ) / (this%timestamps(i) - this%timestamps(i-1))
              exit
            endif
          end do
        endif
      endif
    end select

    ! scale normalized load amplitude by intensity
    real_amplitude = real_amplitude*this%intensity
    
    return
  
  end function evaluateamplitude
  
end module class_amplitude
