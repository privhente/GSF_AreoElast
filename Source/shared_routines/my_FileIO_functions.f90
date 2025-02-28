module my_FileIO_functions
  use SingPrec
  implicit none
  
  contains

!> This routine opens an unformatted output file.
  subroutine OpenBOutFileToRead(Un, OutFile, io_error)
    integer, intent(in)              :: Un
    character(*), intent(in)         :: OutFile                                     !< Name of the output file.
    integer, intent(inout)           :: io_error                                    !< I/O status of OPEN
    character(*), parameter          :: RoutineName = 'OpenFOutFile'

     ! Open output file.  Make sure it worked.
    open(Un, FILE=TRIM(OutFile), STATUS='OLD', FORM='UNFORMATTED' , ACCESS='STREAM', IOSTAT=io_error, ACTION='READ' )
    if (io_error /= 0 )  then
      print*, 'OpenFOutFile:Cannot open file "'//TRIM(OutFile) // '". Another program like MS Excel may have locked it for writing.'
    end if
  
    RETURN
  end subroutine OpenBOutFileToRead  
   
!> This routine opens a formatted output file.
  subroutine OpenFOutFileToRead(Un, OutFile, io_error)
    integer, intent(in)             :: Un
    character(*), intent(in)        :: OutFile                                     !< Name of the output file.
    integer, intent(inout)          :: io_error                                   !< I/O status of OPEN
    character(*), parameter         :: RoutineName = 'OpenFOutFile'

     ! Open output file.  Make sure it worked.
    open( Un, FILE=trim(OutFile), STATUS='replace', FORM='FORMATTED', IOSTAT=io_error, ACTION="write" )
    if (io_error /= 0 )  then
      print*, 'OpenFOutFile:Cannot open file "'//TRIM( OutFile ) // '". Another program like MS Excel may have locked it for writing.'
    end if
  
    RETURN
  end subroutine OpenFOutFileToRead

!< Subroutine for writing integer vector into file
  subroutine writeIntVectorToFile(Sarr,fileID,strfileName)
    integer, Dimension(:) :: Sarr
    integer :: fileID, i
    character(*) :: strfileName

    open(unit = fileID, file = strfileName, action = 'write', status = 'replace')
    do i=1,size(Sarr)
      write(fileID,*) real(Sarr(i))
    end do
    close(fileID)
  end subroutine writeIntVectorToFile
  
!< Subroutine for writing real vector into file
  subroutine writeRealVectorToFile(Sarr,fileID,strfileName)
      real(kind = 8), Dimension(:) :: Sarr
      integer :: fileID, i
      character(*) :: strfileName

1100 format(1000000E27.18E3)
1000 format(1000000f30.15)
      
      open(unit = fileID, file = strfileName, action = 'write', status = 'replace')
      do i=1,size(Sarr)
        write(fileID,1100) Sarr(i)
      end do
      close(fileID)
  end subroutine writeRealVectorToFile

!< Subroutine for writing real matrix into file
  subroutine writeRealMatrixToFile(Sarr,fileID,strfileName)
	  real(kind = 8), Dimension(:,:) :: Sarr
	  integer :: fileID, i, j
	  character(*) :: strfileName
  
1100 format(1000000E27.18E3)
1000 format(1000000f30.15)
     
	  open(unit = fileID, file = strfileName, action = 'write', status = 'replace')
	  do i=1,size(Sarr,1)
		  do j=1,size(Sarr,2)
		    write(fileID,1100,advance='no') Sarr(i,j)
		  end do
		  write(fileID,*) !''  ! this gives you the line break
	  end do
	  close(fileID)
    return
  end subroutine writeRealMatrixToFile

!> This routine returns the next unit number greater than 9 that is not currently in use.
  subroutine GetNewUnit (UnIn)
    integer, INTENT(OUT)     :: UnIn                                         !< Logical unit for the file.
    integer                  :: Un                                           ! Unit number
    logical                  :: Opened                                       ! Flag indicating whether or not a file is opened.
    integer, parameter       :: StartUnit = 10                               ! Starting unit number to check (numbers less than 10 reserved)
    integer, PARAMETER       :: MaxUnit   = 999                               ! The maximum unit number available (or 10 less than the number of files you want to have open at a time)

    ! Initialize subroutine outputs
    Un = StartUnit

    ! See if unit is connected to an open file. Check the next largest number until it is not opened.
    do
      inquire(UNIT=Un, OPENED=Opened)
      if (Opened .eqv. .FALSE.) exit
      Un = Un + 1
      if (Un > MaxUnit) then
        print*,'GetNewUnit() was unable to find an open file unit specifier'
        exit           ! stop searching now
      end if
    end do
    UnIn = Un
    return
   end subroutine GetNewUnit
   
  end module my_FileIO_functions