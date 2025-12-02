!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Generate a soil moisture mask
!>
module generate_moist_mask_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_WRITE,         &
                                      ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2, &
                                      CELL_COLUMN
  use constants_mod,           only : i_def, r_def
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: generate_moist_mask_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! moist_mask
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! unfrozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_sat
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_wilt
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1)  & ! land_fraction
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: generate_moist_mask_code
  end type

  public :: generate_moist_mask_code

contains

  !> @brief Generate a soil moisture mask
  !> @details The lake water conservation correction only applies to grid cells
  !> where the unfrozen soil moisture at the lowest soil level exceeds the
  !> wilting point. This kernel finds out where those points are and returns a
  !> real mask where:
  !> 1.0 = Soil is moist enough to apply lake water correction
  !> 0.0 = Soil is too dry to apply lake water correction
  !> @param[in]     nlayers                Number of layers
  !> @param[in,out] moist_mask             The moisture mask to return
  !> @param[in]     unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in]     soil_moist_sat         Volumetric soil moist at saturation
  !> @param[in]     soil_moist_wilt        Volumetric soil moist at wilting pt
  !> @param[in]     land_fraction          The land fraction
  subroutine generate_moist_mask_code(     &
               nlayers, moist_mask,        &
               unfrozen_soil_moisture,     &
               soil_moist_sat,             &
               soil_moist_wilt,            &
               land_fraction,              &
               ndf_2d,                     &
               undf_2d,                    &
               map_2d,                     &
               ndf_soil,                   &
               undf_soil,                  &
               map_soil                    &
                          )

    use nlsizes_namelist_mod, only: sm_levels

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d, ndf_soil, undf_soil
    integer(kind=i_def), dimension(ndf_2d),    intent(in) :: map_2d
    integer(kind=i_def), dimension(ndf_soil),  intent(in) :: map_soil

    real(kind=r_def), intent(inout) :: moist_mask(undf_2d)
    real(kind=r_def), intent(in)    :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(in)    :: soil_moist_sat(undf_2d)
    real(kind=r_def), intent(in)    :: soil_moist_wilt(undf_2d)
    real(kind=r_def), intent(in)    :: land_fraction(undf_2d)

    real(kind=r_def), parameter     :: tiny = 1.0e-20_r_def

    if ( land_fraction(map_2d(1)) > tiny               .and.                             &
         soil_moist_wilt(map_2d(1)) > tiny             .and.                             &
         unfrozen_soil_moisture(map_soil(1) + sm_levels - 1) * soil_moist_sat(map_2d(1)) >  soil_moist_wilt(map_2d(1)) ) then
      moist_mask(map_2d(1)) = 1.0_r_def
    else
      moist_mask(map_2d(1)) = 0.0_r_def
    end if

  end subroutine generate_moist_mask_code

end module generate_moist_mask_kernel_mod
