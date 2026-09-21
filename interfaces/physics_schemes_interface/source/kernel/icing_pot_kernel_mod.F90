!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Combine relative humidity, temperature and cloud fraction into an
!!        icing-potential index on pressure levels.

module icing_pot_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_SCALAR,       &
                           GH_READ, GH_READWRITE,     &
                           GH_INTEGER,                &
                           GH_REAL, CELL_COLUMN,      &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use conversions_mod, only: zerodegc
  use kernel_mod,    only: kernel_type

  implicit none

  private

  ! Freezing point (0 degrees Celsius) and the lower bound (-20 degrees
  ! Celsius) of the temperature window in which supercooled cloud, and hence
  ! aircraft icing, can occur.
  real(kind=r_def), parameter :: freezing_point = real(zerodegc, r_def)
  real(kind=r_def), parameter :: icing_temp_range = 20.0_r_def
  ! Upper bound used to reject unphysical cloud-fraction values.
  real(kind=r_def), parameter :: max_cloud_fraction = 1.1_r_def
  ! Relative humidity is expressed as a 0-1 fraction and capped at saturation.
  real(kind=r_def), parameter :: max_relative_humidity = 1.0_r_def

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: icing_pot_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR,GH_INTEGER, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: icing_pot_code
  end type icing_pot_kernel_type

  public :: icing_pot_code

contains

  !> @details Forms the icing-potential index as the relative humidity
  !!          (capped at saturation) retained only where cloud is present and
  !!          the temperature lies between 0 and -20 degrees Celsius; it is
  !!          zero elsewhere.
  !> @param[in]     nlayers      The number of layers
  !> @param[in,out] icing_pot    Icing-potential index (dimensionless)
  !> @param[in]     rh           Relative humidity fraction on pressure levels
  !> @param[in]     temperature  Temperature on pressure levels (K)
  !> @param[in]     cloud_frac   Bulk cloud fraction on pressure levels
  !> @param[in]     nplev        Number of pressure levels
  !> @param[in]     ndf          Number of degrees of freedom per cell
  !> @param[in]     undf         Number of total degrees of freedom
  !> @param[in]     map          Dofmap for the cell at the base of the column
  subroutine icing_pot_code(nlayers,     &
                            icing_pot,   &
                            rh,          &
                            temperature, &
                            cloud_frac,  &
                            nplev,       &
                            ndf, undf, map)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nplev
    integer(kind=i_def), intent(in) :: ndf, undf
    integer(kind=i_def), intent(in), dimension(ndf)  :: map

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in),    dimension(undf) :: rh, temperature, &
                                                        cloud_frac
    real(kind=r_def), intent(inout), dimension(undf) :: icing_pot

    ! Internal variables
    integer(kind=i_def) :: kp, idx
    real(kind=r_def)    :: rh_capped, temp_mask, cloud_mask

    do kp = 1, nplev

      idx = map(1) + kp - 1

      ! Cap the relative humidity at saturation.
      rh_capped = min(rh(idx), max_relative_humidity)

      ! Retain only gridboxes where the temperature is between 0 and -20 C.
      temp_mask = merge(1.0_r_def, 0.0_r_def,                             &
                        temperature(idx) > (freezing_point                &
                                            - icing_temp_range)           &
                        .and. temperature(idx) < freezing_point)

      ! Retain only gridboxes where cloud is present.
      cloud_mask = merge(1.0_r_def, 0.0_r_def,                            &
                         cloud_frac(idx) > 0.0_r_def                      &
                         .and. cloud_frac(idx) < max_cloud_fraction)

      icing_pot(idx) = rh_capped * temp_mask * cloud_mask

    end do

  end subroutine icing_pot_code

end module icing_pot_kernel_mod
