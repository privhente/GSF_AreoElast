module IWinAPI_m
! Purpose    : Interfaces toward Microsoft Windows API functions
! Author     : FortranFan
! Reference  : Microsoft Documentation
!
! Description:
! This module defines the interfaces and type aliases toward Windows APIs.
!
   use, intrinsic :: iso_c_binding, only : c_char, c_int, c_intptr_t, c_funptr

   implicit none

   !.. Public by default
   public

   !.. Mnemonics for types in Windows API functions
   integer(c_int), parameter :: HMODULE = c_intptr_t  !.. A handle to a module; base address in memory
   integer(HMODULE), parameter :: NULL_HANDLE = int( 0, kind=kind(NULL_HANDLE) )
   integer(c_int), parameter :: BOOL = c_int

   interface

      function LoadLibrary(lpLibName) bind(C, name='LoadLibraryA') result(RetVal)
      !DIR$ ATTRIBUTES STDCALL :: LoadLibrary
      ! https://docs.microsoft.com/en-us/windows/win32/api/libloaderapi/nf-libloaderapi-loadlibrarya
      ! HMODULE LoadLibraryA( [in] LPCSTR lpLibFileName );

         import :: c_char, HMODULE

         !.. Argument list
         character(kind=c_char, len=1), intent(in) :: lpLibName(*)
         !.. Function result
         integer(HMODULE) :: RetVal

      end function LoadLibrary

      function FreeLibrary(lpLibHandle) bind(C, NAME='FreeLibrary') result(RetVal)
      !DIR$ ATTRIBUTES STDCALL :: FreeLibrary
      ! https://docs.microsoft.com/en-us/windows/win32/api/libloaderapi/nf-libloaderapi-freelibrary
      ! BOOL FreeLibrary( [in] HMODULE hLibModule );

         import :: HMODULE, BOOL

         !.. Argument list
         integer(HMODULE), value :: lpLibHandle
         !.. Function result
         integer(BOOL) :: RetVal

      end function FreeLibrary

      function GetProcAddress(lpLibHandle, lpProcName) bind(C, name='GetProcAddress') result(RetVal)
      !DIR$ ATTRIBUTES STDCALL :: GetProcAddress
      ! https://docs.microsoft.com/en-us/windows/win32/api/libloaderapi/nf-libloaderapi-getprocaddress
      ! FARPROC GetProcAddress( [in] HMODULE hModule, [in] LPCSTR  lpProcName );

         import :: HMODULE, c_char, c_funptr

         !.. Argument list
         integer(HMODULE), intent(in), value      :: lpLibHandle
         character(kind=c_char,len=1), intent(in) :: lpProcName(*)
         !.. Function result
         type(c_funptr) :: RetVal

      end function GetProcAddress

   end interface
    
end module IWinAPI_m
