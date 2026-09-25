!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Icing potential on pressure levels.

module icing_pot_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_SCALAR,       &
                           GH_READ, GH_WRITE,         &
                           GH_REAL, GH_INTEGER,       &
                           CELL_COLUMN,               &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: icing_pot_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                          &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: icing_pot_code
  end type icing_pot_kernel_type

  public :: icing_pot_code

contains

  !> @brief   Calculates the icing potential on pressure levels.
  !> @details The icing potential is the relative humidity, capped at one,
  !>          where cloud is present and the temperature is between -20 and
  !>          0 degrees Celsius, and zero elsewhere. All inputs are already
  !>          on pressure levels.
  !> @param[in]     nlayers     Number of layers
  !> @param[in,out] icing_pot   Icing potential on pressure levels (fraction)
  !> @param[in]     rh          Relative humidity on pressure levels
  !> @param[in]     temperature Temperature on pressure levels (K)
  !> @param[in]     cloud       Bulk cloud fraction on pressure levels
  !> @param[in]     nplev       Number of pressure levels
  !> @param[in]     zerodegc    Zero degrees Celsius in Kelvin
  !> @param[in]     ndf         Number of degrees of freedom per cell
  !> @param[in]     undf        Number of total degrees of freedom
  !> @param[in]     map         Dofmap for the cell at the base of the column
  subroutine icing_pot_code(nlayers,        &
                            icing_pot,      &
                            rh,             &
                            temperature,    &
                            cloud,          &
                            nplev,          &
                            zerodegc,       &
                            ndf, undf, map)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nplev
    integer(kind=i_def), intent(in) :: ndf, undf
    integer(kind=i_def), intent(in), dimension(ndf) :: map

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(inout), dimension(undf) :: icing_pot
    real(kind=r_def), intent(in),    dimension(undf) :: rh, temperature, &
                                                        cloud
    real(kind=r_def), intent(in) :: zerodegc

    ! Width of the icing temperature range below zero degrees Celsius (K)
    real(kind=r_def), parameter :: icing_temp_range = 20.0_r_def
    ! Maximum relative humidity used in the icing potential
    real(kind=r_def), parameter :: rh_max = 1.0_r_def
    ! Upper bound of the valid cloud fraction range
    real(kind=r_def), parameter :: cloud_max = 1.1_r_def

    ! Internal variables
    integer(kind=i_def) :: kp, df
    real(kind=r_def)    :: temp_mask, cloud_mask

    do kp = 1, nplev
      df = map(1) + kp - 1

      ! Temperature between -20 and 0 degrees Celsius
      temp_mask = merge( 1.0_r_def, 0.0_r_def,                         &
                         temperature(df) > zerodegc - icing_temp_range &
                         .and. temperature(df) < zerodegc )

      ! Cloud present
      cloud_mask = merge( 1.0_r_def, 0.0_r_def,                        &
                          cloud(df) > 0.0_r_def .and. cloud(df) < cloud_max )

      icing_pot(df) = min( rh(df), rh_max ) * cloud_mask * temp_mask
    end do

  end subroutine icing_pot_code

end module icing_pot_kernel_mod
