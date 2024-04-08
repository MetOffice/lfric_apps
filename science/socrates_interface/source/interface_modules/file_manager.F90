!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief File management routines for assignment and release of file units
! @detail Service routines provided to Socrates by LFRic (io_utility_mod)

module file_manager

  use constants_mod,  only: i_def
  use io_utility_mod, only: claim_io_unit, release_io_unit

  implicit none
  private
  public :: assign_file_unit, release_file_unit

contains

  subroutine assign_file_unit(filename, iunit, handler)

    implicit none
    character(*),   intent(in)  :: filename
    integer(i_def), intent(out) :: iunit
    character(*),   intent(in)  :: handler

    iunit = claim_io_unit()

  end subroutine assign_file_unit

  subroutine release_file_unit(iunit, handler)

    implicit none
    integer(i_def), intent(in) :: iunit
    character(*),   intent(in) :: handler

    integer(i_def) :: dummy_iunit

    call release_io_unit(dummy_iunit)

  end subroutine release_file_unit

end module file_manager
