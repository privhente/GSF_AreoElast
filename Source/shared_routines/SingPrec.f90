!< module for setting parameter for defined type
module SingPrec
  implicit none
  integer, parameter  :: B2Ki = SELECTED_INT_KIND(4)         !< Kind for two-byte whole numbers
  integer, parameter  :: B4Ki = SELECTED_INT_KIND(9)         !< Kind for four-byte whole number
  integer, parameter  :: SiKi = SELECTED_REAL_KIND(6,30)     !< Kind for four-byte, floating-point numbers
end module SingPrec