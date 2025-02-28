module DllHelper_m
! Purpose    : Helper module and type to work with Dlls on Windows
! Author     : FortranFan
! Reference  : Using Run-Time Dynamic Linking
!              https://docs.microsoft.com/en-us/windows/win32/dlls/using-run-time-dynamic-linking
! Description:
! This module defines a public type, DllHelper_t, to work with DLLs on Windows.  This type has
!  - Load method with an Dll name as input to load the DLL
!  - GetFunPtr to dispatch a C function pointer for an exported procedure in the Dll; the proc
!    name string is the input
!  - a finalizer that frees up the Dll when the object is destroyed.
!
   use, intrinsic :: iso_c_binding, only : c_char, c_funptr, c_null_funptr
   use IWINApi_m, only : HMODULE, NULL_HANDLE, BOOL, LoadLibrary, GetProcAddress, FreeLibrary

   private

   type, public :: DllHelper_t
      private
      character(kind=c_char, len=:), allocatable :: m_DllName
      integer(HMODULE) :: m_DllHandle = NULL_HANDLE
   contains
      final :: FreeDll
      procedure, pass(this) :: Load => LoadDll
      procedure, pass(this) :: GetFunPtr
   end type
   
contains
    
   subroutine FreeDll( this )
   ! Unload the DLL and free up its resources
      type(DllHelper_t), intent(inout) :: this
      
      ! Local variables
      integer(BOOL) :: iret
      
      if ( this%m_DllHandle /= NULL_HANDLE ) then
         iret = FreeLibrary( this%m_DllHandle )
      end if
      this%m_DllHandle = NULL_HANDLE
      
      return
        
   end subroutine 

   subroutine LoadDll( this, DllName, iret )
   ! Load the DLL and set up its handle

      class(DllHelper_t), intent(inout) :: this
      character(kind=c_char, len=*), intent(in) :: DllName
      integer(BOOL), intent(inout) :: iret

      iret = 0
      this%m_DllHandle = LoadLibrary( DllName )
      if ( this%m_DllHandle == NULL_HANDLE ) then
         iret = 1 !<-- TODO: replace with GetLastError
      end if
      
      return
        
   end subroutine 

   function GetFunPtr( this, ProcName ) result(FunPtr)
   ! Get the function pointer

      class(DllHelper_t), intent(inout) :: this
      character(kind=c_char, len=*), intent(in) :: ProcName
      ! Function result
      type(c_funptr) :: FunPtr

      FunPtr = GetProcAddress( this%m_DllHandle, ProcName )
      
      return
        
   end function 

end module 