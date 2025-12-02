!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Converts the vertical departure distance into a fractional wind flux.
!> @details The Eulerian departure distance is calculated using the advecting
!!          wind divided by Det(J) multiplied by the time step. This kernels
!!          splits the vertical departure distance into integer and fractional
!!          parts, and then returns the fractional vertical advecting wind flux,
!!          which will have units of volume.

module fractional_vertical_wind_kernel_mod

use argument_mod,       only : arg_type,              &
                               GH_FIELD, GH_REAL,     &
                               GH_READ, GH_WRITE,     &
                               CELL_COLUMN
use fs_continuity_mod,  only : W2v, W3
use constants_mod,      only : r_tran, i_def, EPS_R_TRAN
use kernel_mod,         only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: fractional_vertical_wind_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! frac_wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep_dist
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3)   & ! detj
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: fractional_vertical_wind_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: fractional_vertical_wind_code

contains

  !> @brief Compute the flux using the PCM reconstruction.
  !> @param[in]     nlayers   Number of layers
  !> @param[in,out] frac_wind The fractional vertical wind
  !> @param[in]     dep_dist  The vertical departure distance
  !> @param[in]     detj      Det(J) at W3 points
  !> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
  !> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
  !> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
  !> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3   Number of unique degrees of freedom for W3
  !> @param[in]     map_w3    Map for W3
  subroutine fractional_vertical_wind_code( nlayers,   &
                                            frac_wind, &
                                            dep_dist,  &
                                            detj,      &
                                            ndf_w2v,   &
                                            undf_w2v,  &
                                            map_w2v,   &
                                            ndf_w3,    &
                                            undf_w3,   &
                                            map_w3 )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: undf_w2v
    integer(kind=i_def), intent(in)    :: ndf_w2v
    integer(kind=i_def), intent(in)    :: undf_w3
    integer(kind=i_def), intent(in)    :: ndf_w3
    real(kind=r_tran),   intent(inout) :: frac_wind(undf_w2v)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
    real(kind=r_tran),   intent(in)    :: detj(undf_w3)
    integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
    integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

    ! Internal variables
    integer(kind=i_def) :: k, df, km, kp

    real(kind=r_tran)   :: departure_dist
    real(kind=r_tran)   :: fractional_distance
    real(kind=r_tran)   :: upwind_detj
    integer(kind=i_def) :: int_dist

    ! Wind at surface and model top should be zero
    frac_wind(map_w2v(1)) = 0.0_r_tran
    frac_wind(map_w2v(2)+nlayers-1) = 0.0_r_tran

    ! Loop over levels working with DOF at the top of the cell
    df = 2

    do k = 0, nlayers-2
      ! Get departure distance
      departure_dist = dep_dist(map_w2v(df)+k)

      ! Get integer and fractional part
      fractional_distance = departure_dist - int(departure_dist)
      int_dist = int(departure_dist)

      ! Get indices
      km = max( 0_i_def, min(k - int_dist, nlayers-1) )
      kp = max( 0_i_def, min(k + 1 - int_dist, nlayers-1) )

      ! Compute upwind Det(J)
      upwind_detj = ( 0.5_r_tran + sign(0.5_r_tran, departure_dist) ) * detj(map_w3(1) + km) + &
                    ( 0.5_r_tran - sign(0.5_r_tran, departure_dist) ) * detj(map_w3(1) + kp)

      ! Fractional wind is fractional_departure_distance * Det(J)
      frac_wind(map_w2v(df)+k) = fractional_distance * upwind_detj
    end do

  end subroutine fractional_vertical_wind_code

end module fractional_vertical_wind_kernel_mod
