!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Applies the lake water correction to soil moisture
!>
module apply_lake_correction_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_WRITE,         &
                                      GH_SCALAR,                 &
                                      ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2, &
                                      CELL_COLUMN
  use constants_mod,           only : i_def, r_def
  use kernel_mod,              only : kernel_type
  use driver_water_constants_mod, only: density_h2o

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: apply_lake_correction_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! moist_mask
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! unfrozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_sat
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! soil_moisture
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                 & ! water_reduction
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: apply_lake_correction_code
  end type

  public :: apply_lake_correction_code

contains

  !> @brief Applies the lake water correction to soil moisture
  !> @details Reduce the soil moisture at moist points (also at lowest soil level)
  !> to account for lake evaporation
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     moist_mask             The mask where soil moisture is moist emough to reduce it
  !> @param[in,out] unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in]     soil_moist_sat         Volumetric soil moist at saturation
  !> @param[in,out] soil_moisture          Soil moisture content (kg m-2)
  !> @param[in]     water_reduction        The amount of water to reduce (kg m-2)
  subroutine apply_lake_correction_code(   &
               nlayers, moist_mask,        &
               unfrozen_soil_moisture,     &
               soil_moist_sat,             &
               soil_moisture,              &
               water_reduction,            &
               ndf_2d,                     &
               undf_2d,                    &
               map_2d,                     &
               ndf_soil,                   &
               undf_soil,                  &
               map_soil                    &
                          )

    use nlsizes_namelist_mod, only: sm_levels
    use jules_soil_mod, only: dzsoil

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d, ndf_soil, undf_soil
    integer(kind=i_def), dimension(ndf_2d),    intent(in) :: map_2d
    integer(kind=i_def), dimension(ndf_soil),  intent(in) :: map_soil

    real(kind=r_def), intent(in)    :: moist_mask(undf_2d)
    real(kind=r_def), intent(inout) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(in)    :: soil_moist_sat(undf_2d)
    real(kind=r_def), intent(inout) :: soil_moisture(undf_soil)
    real(kind=r_def), intent(in)    :: water_reduction

    if ( moist_mask(map_2d(1)) >  0.5_r_def ) then
      soil_moisture(map_soil(1) + sm_levels - 1 ) = soil_moisture(map_soil(1) + sm_levels - 1 ) - water_reduction
      unfrozen_soil_moisture(map_soil(1) + sm_levels - 1 ) = unfrozen_soil_moisture(map_soil(1) + sm_levels - 1 ) &
                                   - water_reduction / (soil_moist_sat(map_2d(1)) * density_h2o * dzsoil(sm_levels))
    end if

  end subroutine apply_lake_correction_code

end module apply_lake_correction_kernel_mod
