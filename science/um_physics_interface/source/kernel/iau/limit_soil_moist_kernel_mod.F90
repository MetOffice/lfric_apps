!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Constrain updated soil moisture
module limit_soil_moist_kernel_mod

  use argument_mod,                only: arg_type,                  &
                                         GH_FIELD, GH_SCALAR,       &
                                         GH_REAL, GH_INTEGER,       &
                                         GH_READ, GH_READWRITE,     &
                                         ANY_DISCONTINUOUS_SPACE_1, &
                                         ANY_DISCONTINUOUS_SPACE_2, &
                                         CELL_COLUMN
  use constants_mod,               only: i_def, r_def, RMDI, i_um, r_um
  use kernel_mod,                  only: kernel_type
  use jules_soil_mod,              only: dzsoil

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: limit_soil_moist_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),      &
         arg_type(GH_FIELD,  GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),      &
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: limit_soil_moist_code
  end type limit_soil_moist_kernel_type

  public :: limit_soil_moist_code

contains

  !> @param[in]     nlayers          Number of layers
  !> @param[in]     soil_moist_wilt  Volumetric soil moist at wilting pt
  !> @param[in]     soil_moist_sat   Volumetric soil moist at saturation
  !> @param[in,out] soil_moisture    Soil moisture on soil levels
  !> @param[in]     sm_levels        Number of soil moisture levels
  !> @param[in]     rho_water        Water density
  !> @param[in]     ndf_2d           Number of total DOFs per cell for 2d fields
  !> @param[in]     undf_2d          Number of total DOFs for 2d fields
  !> @param[in]     map_2d           Dofmap for cell for surface 2d fields
  !> @param[in]     ndf_soil         Number of DOFs per call for soil levels
  !> @param[in]     undf_soil        Number of total DOFs for soil levels
  !> @param[in]     map_soil         Dofmap for cell for soil levels
  subroutine limit_soil_moist_code ( nlayers,                      &
                                     soil_moist_wilt,              &
                                     soil_moist_sat,               &
                                     soil_moisture,                &
                                     sm_levels,                    &
                                     rho_water,                    &
                                     ndf_2d, undf_2d, map_2d,      &
                                     ndf_soil, undf_soil, map_soil )

    implicit none

    ! Arguments
    integer( kind=i_def ), intent(in)    :: nlayers
    integer( kind=i_def ), intent(in)    :: sm_levels
    real( kind=r_def ),    intent(in)    :: rho_water

    integer( kind=i_def ), intent(in)    :: ndf_2d, undf_2d
    integer( kind=i_def ), intent(in)    :: map_2d(ndf_2d)
    integer( kind=i_def ), intent(in)    :: ndf_soil, undf_soil
    integer( kind=i_def ), intent(in)    :: map_soil(ndf_soil)

    real( kind=r_def ),    intent(in)    :: soil_moist_wilt(undf_2d)
    real( kind=r_def ),    intent(in)    :: soil_moist_sat(undf_2d)
    real( kind=r_def ),    intent(inout) :: soil_moisture(undf_soil)

    ! local variables
    integer( i_um ) :: n
    real( r_um ) :: zdepth, smc_min, smc_max

    ! check soil moisture against 0.1*wilting point and saturation point limits
    ! and constrain to these limits

    ! only process valid soil points
    if ( soil_moist_sat(map_2d(1)) > 0.0_r_def ) then

      do n = 1, sm_levels

        ! calculate the limits for each soil level
        smc_min = 0.1_r_def * soil_moist_wilt(map_2d(1)) * dzsoil(n) * &
                  rho_water

        smc_max = soil_moist_sat(map_2d(1)) * dzsoil(n) * rho_water

        ! now apply the caulculated limits
        if ( soil_moisture(map_soil(1) + n-1) < smc_min ) then
          soil_moisture(map_soil(1) + n-1) = smc_min
        else if ( soil_moisture(map_soil(1) + n-1) > smc_max ) then
          soil_moisture(map_soil(1) + n-1) = smc_max
        end if

      end do

    end if

  end subroutine limit_soil_moist_code

end module limit_soil_moist_kernel_mod
