module DISCON_m
! Module to define Fortran interfaces of procedures to be consumed from C++ Dll(s) 

   use, intrinsic :: iso_c_binding, only : c_float, c_int, c_char
   use DllHelper_m, only : DllHelper_t

   private
   
   abstract interface 
      subroutine IDISCON ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG ) bind(C)
      ! Assumes the C++ function prototype is as follows:
      ! extern "C" void DISCON(float *, int *, const char *, const char *, char *)
         import :: c_float, c_int, c_char
         real(c_float),          intent(inout) :: avrSWAP(*)      
         integer(c_int),         intent(inout) :: aviFAIL        
         character(kind=c_char), intent(in)    :: accINFILE(*)      
         character(kind=c_char), intent(in)    :: avcOUTNAME(*)     
         character(kind=c_char), intent(inout) :: avcMSG(*)        
      end subroutine IDISCON  
   end interface
   
   ! Program variables
   type(DllHelper_t), save :: controller
   procedure(IDISCON), pointer, save, protected, public :: DISCON => null()
   
   public :: SetupDISCON
   
contains

   subroutine SetupDISCON(boolean_abort)
      use, intrinsic :: iso_c_binding, only : c_funptr, c_associated, c_f_procpointer
      
      implicit none
      
      type(c_funptr) :: pDISCON
      integer :: irc
      logical, intent(inout) :: boolean_abort
      
      call Controller%Load( "libdiscon", irc )  !    Changed to the name of the actual DLL
      if ( irc /= 0 ) then
        print*, 'Warning: failed to load libdiscon.dll...'
        boolean_abort = .TRUE.
        return
      end if
      
      pDISCON = Controller%GetFunPtr( "DISCON" )
      if ( .not. c_associated(pDISCON) ) then
         print*, 'Warning: failed to get the address of DISCON procedure'
         boolean_abort = .TRUE.
         return
      end if
      call c_f_procpointer( cptr=pDISCON, fptr=DISCON )
      
      if ( .not. associated(DISCON) ) then
         print*, 'Warning: DISCON is not associated.'
         boolean_abort = .TRUE.
         return
      end if
      
   end subroutine
   
end module