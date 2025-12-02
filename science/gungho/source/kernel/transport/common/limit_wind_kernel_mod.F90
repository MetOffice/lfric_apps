!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Limits wind velocity.
!>
module limit_wind_kernel_mod

  use argument_mod,            only : arg_type,              &
                                      GH_FIELD, GH_REAL,     &
                                      GH_READ, GH_READWRITE, &
                                      GH_SCALAR, DOF
  use constants_mod,           only : r_def
  use fs_continuity_mod,       only : W2
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: limit_wind_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, W2), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      W2), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)           &
         /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: limit_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: limit_wind_code

contains

!> @brief Limits the wind dofs by some measure of CFL limit
!! @param[in,out] wind Wind
!! @param[in] dJ_on_w2 detJ evaluated on w2 points
!! @param[in] dt The model timestep length
!! @param[in] max_cfl Maximum advective Courant number
subroutine limit_wind_code(wind, dJ_on_w2, dt, max_cfl)

  implicit none

  ! Arguments
  real(kind=r_def), intent(inout) :: wind
  real(kind=r_def), intent(in)    :: dJ_on_w2
  real(kind=r_def), intent(in)    :: max_cfl
  real(kind=r_def), intent(in)    :: dt

  ! Internal variables
  real(kind=r_def)    :: wind_max

  ! Set maximum wind flux (i.e. wind speed * face area) to
  ! dJ*maxcfl*dt, i.e. wind speed < dJ/dA *maxcfl/dt ~ dx * maxcfl/dt
  ! dJ~cell volume
  wind_max = dJ_on_w2*max_cfl/dt

  ! Limit both positive and negative directions
  wind = max(min(wind, wind_max), -wind_max)

end subroutine limit_wind_code

end module limit_wind_kernel_mod
