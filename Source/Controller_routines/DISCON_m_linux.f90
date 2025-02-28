module DISCON_m_linux

  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_char

  !< Declare the interface to the subroutine
  interface
    subroutine DISCON( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG ) bind(c, name="DISCON")
       import :: c_float, c_int, c_char
       real(c_float),          intent(inout) :: avrSWAP(*)      
       integer(c_int),         intent(inout) :: aviFAIL        
       character(kind=c_char), intent(in)    :: accINFILE(*)      
       character(kind=c_char), intent(in)    :: avcOUTNAME(*)     
       character(kind=c_char), intent(inout) :: avcMSG(*)  
    end subroutine DISCON
  end interface

end module DISCON_m_linux
