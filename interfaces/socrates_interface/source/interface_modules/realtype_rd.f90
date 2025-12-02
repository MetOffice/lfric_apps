!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Set the precision of real variables for Socrates

module realtype_rd

  use constants_mod, only: r_def, r_um

  implicit none
  private

  ! Internal Socrates precision
  integer, public, parameter :: RealK = r_um

  ! Precision of variables passed through the Socrates interface
  integer, public, parameter :: RealExt = r_def

end module realtype_rd
