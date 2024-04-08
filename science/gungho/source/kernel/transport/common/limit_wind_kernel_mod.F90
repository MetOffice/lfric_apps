!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Limits wind velocity.
!>
module limit_wind_kernel_mod

  use argument_mod,            only : arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_READ, GH_INC,   &
                                      GH_SCALAR, CELL_COLUMN
  use constants_mod,           only : i_def, r_def
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
    type(arg_type) :: meta_args(4) = (/             &
         arg_type(GH_FIELD,  GH_REAL, GH_INC,  W2), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, W2), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),     &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)      &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: limit_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: limit_wind_code

contains

!> @brief Limits the wind dofs by some measure of CFL limit
!! @param[in] nlayers Number of layers
!! @param[in,out] wind Wind
!! @param[in] dJ_on_w2 detJ evaluated on w2 points
!! @param[in] dt The model timestep length
!! @param[in] max_cfl Maximum advective Courant number
!! @param[in] ndf_w2 Number of degrees of freedom per cell
!! @param[in] undf_w2 Total number of degrees of freedom
!! @param[in] map_w2 Dofmap for the cell at the base of the column
subroutine limit_wind_code(nlayers, wind, dJ_on_w2, dt, &
                           max_cfl, ndf_w2, undf_w2, map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf_w2, undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)  :: map_w2
  real(kind=r_def), dimension(undf_w2), intent(inout) :: wind
  real(kind=r_def), dimension(undf_w2), intent(in)    :: dJ_on_w2
  real(kind=r_def),                     intent(in)    :: max_cfl
  real(kind=r_def),                     intent(in)    :: dt

  ! Internal variables
  real(kind=r_def)    :: wind_max
  integer(kind=i_def) :: df, k

  do k = 0, nlayers-1
    do df = 1, ndf_w2

      ! Set maximum wind flux (i.e. wind speed * face area) to
      ! dJ*maxcfl*dt, i.e. wind speed < dJ/dA *maxcfl/dt ~ dx * maxcfl/dt
      ! dJ~cell volume
      wind_max = dJ_on_w2(map_w2(df)+k)*max_cfl/dt

      ! Limit both positive and negative directions
      wind(map_w2(df)+k) = max(min(wind(map_w2(df)+k), wind_max), -wind_max)

    end do
  end do

end subroutine limit_wind_code

end module limit_wind_kernel_mod
