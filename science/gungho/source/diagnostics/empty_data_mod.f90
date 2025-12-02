!-------------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Provides data fields that can be used by local diagnostics
!> @details If a local diagnostic is not requested for output, the field that holds
!>          it can be overridden by passing in an "empty" array, which will mean
!>          that the field will use less memory.
!>          Within the kernel, the data passed in should be checked: if it is
!>          associated with this empty array, the field should not be written to.

module empty_data_mod
  USE constants_mod, ONLY: r_def, i_def

  implicit none
  private
  real(r_def),    target, public ::empty_real_data(1)
  integer(i_def), target, public ::empty_integer_data(1)
end module empty_data_mod
