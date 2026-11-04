!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Met Office Claude Code Enterprise (Claude Opus 5.5).

!> @brief Find a pressure on a pressure-level axis.
!> @details This belongs with get_axis_values in lfric_core's
!>          lfric_xios_diag_mod, and should move there; it lives here until
!>          lfric_core has it.

module plev_index_mod

  use constants_mod, only: r_def, i_def
  use log_mod,       only: log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private

  public :: plev_index

  !> Tolerance for matching a pressure to an axis value (Pa)
  real(kind=r_def), parameter :: plev_tol = 0.1_r_def

contains

  !> @brief   Slot of a pressure on a pressure-level axis.
  !> @details The first axis value within plev_tol of the pressure. A
  !>          pressure not on the axis is an error.
  !> @param[in] axis_name  Name of the axis, for the error message
  !> @param[in] plevs      Values of the axis (Pa)
  !> @param[in] pressure   Pressure to find (Pa)
  !> @return    Slot of the pressure on the axis
  function plev_index(axis_name, plevs, pressure) result(slot)

    implicit none

    character(len=*), intent(in) :: axis_name
    real(kind=r_def), intent(in) :: plevs(:)
    real(kind=r_def), intent(in) :: pressure
    integer(kind=i_def)          :: slot

    integer(kind=i_def) :: k

    slot = 0_i_def
    do k = 1, size(plevs)
      if ( abs(plevs(k) - pressure) < plev_tol ) then
        slot = k
        exit
      end if
    end do

    if ( slot == 0_i_def ) then
      write(log_scratch_space, '(A, F10.1, A, A)')                          &
           'Pressure level ', pressure, ' Pa is not on axis ',              &
           trim(axis_name)
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end function plev_index

end module plev_index_mod
