!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Perform IAU updates for surface increments
module surf_inc_update_kernel_mod

  use argument_mod,               only: arg_type,                  &
                                        GH_FIELD, GH_SCALAR,       &
                                        GH_REAL, GH_INTEGER,       &
                                        GH_READ, GH_READWRITE,     &
                                        ANY_DISCONTINUOUS_SPACE_1, &
                                        ANY_DISCONTINUOUS_SPACE_2, &
                                        ANY_DISCONTINUOUS_SPACE_3, &
                                        ANY_DISCONTINUOUS_SPACE_4, &
                                        ANY_DISCONTINUOUS_SPACE_5, &
                                        CELL_COLUMN
  use constants_mod,              only: i_def, r_def, l_def, r_um
  use kernel_mod,                 only: kernel_type


  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: surf_inc_update_kernel_type
    private
    type(arg_type) :: meta_args(15) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & !soil_moist_sat
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & !grid_snow_mass
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & !soil_temperature_inc
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & !soil_temperature
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & !soil_moisture_inc
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & !soil_moisture
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & !tile_temperature_inc
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & !tile_temperature
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & !snow_layer_temp_inc
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & !snow_layer_temp
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !n_land_tile
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !snow_lev_tile
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !nsmax
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !sm_levels
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  & !T_freeze_h2o
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: surf_inc_update_code
  end type surf_inc_update_kernel_type

  public :: surf_inc_update_code

contains

  !> @param[in]     nlayers              Number of layers
  !> @param[in]     soil_moist_sat       Volumetric soil moist at saturation
  !> @param[in]     grid_snow_mass       Grid-box mean snow mass
  !> @param[in]     soil_temperature_inc Soil temp increment on soil levels
  !> @param[in,out] soil_temperature     Soil temperature on soil levels
  !> @param[in]     soil_moisture_inc    Soil moisture increment on soil levels
  !> @param[in,out] soil_moisture        Soil moisture on soil levels
  !> @param[in]     tile_temperature_inc Tile temp increment on land tiles
  !> @param[in,out] tile_temperature     Tile temperature on surface tiles
  !> @param[in]     snow_layer_temp_inc  Snow temp increment on snow layers and tiles
  !> @param[in,out] snow_layer_temp      Snow temp on snow layers and tiles
  !> @param[in]     n_land_tile          Number of land tiles
  !> @param[in]     snow_lev_tile        Number of snow levels and land tiles
  !> @param[in]     sm_levels            Number of soil levels
  !> @param[in]     nsmax                Number of multilevels in snow scheme
  !> @param[in]     T_freeze_h2o         Freezing point of water
  !> @param[in]     ndf_2d               number of total DOFs per cell for 2d fields
  !> @param[in]     undf_2d              number of total DOFs for 2d fields
  !> @param[in]     map_2d               dofmap for cell for surface 2d fields
  !> @param[in]     ndf_soil             number of DOFs per call for soil levels
  !> @param[in]     undf_soil            number of total DOFs for soil levels
  !> @param[in]     map_soil             Dofmap for cell for soil levels
  !> @param[in]     ndf_land             number of total DOFs per cell for land fields
  !> @param[in]     undf_land            number of total DOFs for land fields
  !> @param[in]     map_land             dofmap for cell for surface land fields
  !> @param[in]     ndf_tile             number of DOFs per call for surface fields
  !> @param[in]     undf_tile            number of total DOFs for surface fields
  !> @param[in]     map_tile             Dofmap for cell for surface fields
  !> @param[in]     ndf_snow             number of DOFs per call for snow fields
  !> @param[in]     undf_snow            number of total DOFs for snow fields
  !> @param[in]     map_snow             Dofmap for cell for snow fields

  subroutine surf_inc_update_code ( nlayers,                       &
                                    soil_moist_sat,                &
                                    grid_snow_mass,                &
                                    soil_temperature_inc,          &
                                    soil_temperature,              &
                                    soil_moisture_inc,             &
                                    soil_moisture,                 &
                                    tile_temperature_inc,          &
                                    tile_temperature,              &
                                    snow_layer_temp_inc,           &
                                    snow_layer_temp,               &
                                    n_land_tile,                   &
                                    snow_lev_tile,                 &
                                    sm_levels,                     &
                                    nsmax,                         &
                                    T_freeze_h2o,                  &
                                    ndf_2d, undf_2d, map_2d,       &
                                    ndf_soil, undf_soil, map_soil, &
                                    ndf_land, undf_land, map_land, &
                                    ndf_tile, undf_tile, map_tile, &
                                    ndf_snow, undf_snow, map_snow )

    implicit none

    ! Arguments
    integer( kind=i_def ), intent(in)    :: nlayers
    integer( kind=i_def ), intent(in)    :: n_land_tile
    integer( kind=i_def ), intent(in)    :: snow_lev_tile
    integer( kind=i_def ), intent(in)    :: sm_levels
    integer( kind=i_def ), intent(in)    :: nsmax
    real( kind=r_def ),    intent(in)    :: T_freeze_h2o

    integer( kind=i_def ), intent(in)    :: ndf_2d, undf_2d
    integer( kind=i_def ), intent(in)    :: map_2d(ndf_2d)
    integer( kind=i_def ), intent(in)    :: ndf_soil, undf_soil
    integer( kind=i_def ), intent(in)    :: map_soil(ndf_soil)
    integer( kind=i_def ), intent(in)    :: ndf_land, undf_land
    integer( kind=i_def ), intent(in)    :: map_land(ndf_land)
    integer( kind=i_def ), intent(in)    :: ndf_tile, undf_tile
    integer( kind=i_def ), intent(in)    :: map_tile(ndf_tile)
    integer( kind=i_def ), intent(in)    :: ndf_snow, undf_snow
    integer( kind=i_def ), intent(in)    :: map_snow(ndf_snow)

    real( kind=r_def ),    intent(in)    :: soil_moist_sat(undf_2d)
    real( kind=r_def ),    intent(in)    :: grid_snow_mass(undf_2d)
    real( kind=r_def ),    intent(in)    :: soil_temperature_inc(undf_soil)
    real( kind=r_def ),    intent(inout) :: soil_temperature(undf_soil)
    real( kind=r_def ),    intent(in)    :: soil_moisture_inc(undf_soil)
    real( kind=r_def ),    intent(inout) :: soil_moisture(undf_soil)
    real( kind=r_def ),    intent(in)    :: tile_temperature_inc(undf_land)
    real( kind=r_def ),    intent(inout) :: tile_temperature(undf_tile)
    real( kind=r_def ),    intent(in)    :: snow_layer_temp_inc(undf_snow)
    real( kind=r_def ),    intent(inout) :: snow_layer_temp(undf_snow)

    ! Internal variables
    real( kind=r_um ),     parameter     :: um_rmdi = -1073741824.0_r_um
    integer( kind=i_def )                :: n
    logical( kind=l_def )                :: l_update_soil_temperature(sm_levels)
    logical( kind=l_def )                :: l_update_soil_moisture(sm_levels)
    logical( kind=l_def )                :: l_update_tile_temperature(n_land_tile)
    logical( kind=l_def )                :: l_update_snow_layer_temp(snow_lev_tile)


    ! Update each prognostic field in turn, only if neither increment field
    ! or prognostic field is RMDI. UM and LFRic RMDI values are different and
    ! both may be present

    ! First set logical switches for whether increment should be used
    !-------------------------------------------------------------------

    l_update_soil_temperature(:) = .false.
    l_update_soil_moisture(:) = .false.
    l_update_tile_temperature(:) = .false.
    l_update_snow_layer_temp(:) = .false.

    ! Soil temperature inc only used if holds valid data and not land ice
    ! Soil moisture inc only used if holds valid data, not land ice, and snow
    ! depth <= 0.05 kg/m**2
    if ( soil_moist_sat(map_2d(1)) > 0.0_r_def ) then
      do n = 1, sm_levels
        if ( soil_temperature_inc(map_soil(1)+n-1) > um_rmdi ) l_update_soil_temperature(n) = .true.
        if ( ( soil_moisture_inc(map_soil(1)+n-1) > um_rmdi ) .and. &
          ( grid_snow_mass(map_2d(1)) <= 0.05_r_def) ) then
          l_update_soil_moisture(n) = .true.
        end if
      end do
    end if

    ! Tile temperature inc only used for valid data points on land tiles
    do n = 1, n_land_tile
      if ( ( tile_temperature_inc(map_land(1)+n-1) > um_rmdi ) .and. &
        ( tile_temperature(map_tile(1)+n-1) > um_rmdi ) ) then
        l_update_tile_temperature(n) = .true.
      end if
    end do

    ! Snow layer temperature inc only used for valid data points and when
    ! multilayer snow scheme is active
    if (nsmax > 0) then
      do n = 1, snow_lev_tile
        if ( ( snow_layer_temp_inc(map_snow(1)+n-1) > um_rmdi ) .and. &
          ( snow_layer_temp(map_snow(1)+n-1) > um_rmdi ) ) then
          l_update_snow_layer_temp(n) = .true.
        end if
      end do
    end if

    ! Update each variable in turn
    !---------------------------------

    ! Update variables on soil levels - soil moisture and soil temperature
    do n = 1, sm_levels

      if (l_update_soil_temperature(n)) then
        soil_temperature(map_soil(1)+n-1) = soil_temperature(map_soil(1)+n-1) + &
        soil_temperature_inc(map_soil(1)+n-1)
      end if
      if (l_update_soil_moisture(n)) then
        soil_moisture(map_soil(1)+n-1) = soil_moisture(map_soil(1)+n-1) + &
        soil_moisture_inc(map_soil(1)+n-1)
      end if

    end do

    ! Update variables on surface tiles - tile temperature
    ! Only update land tiles
    do n = 1, n_land_tile

      if (l_update_tile_temperature(n)) then
        tile_temperature(map_tile(1)+n-1) = tile_temperature(map_tile(1)+n-1) + &
        tile_temperature_inc(map_land(1)+n-1)
      end if

    end do

    ! Update variables on snow layers and land tiles - snow layer temperature
    ! Only update if multi-layer snow scheme is active and constrain
    ! to a max value of melting point
    if (nsmax > 0) then

      do n = 1, snow_lev_tile

        if (l_update_snow_layer_temp(n)) then
          snow_layer_temp(map_snow(1)+n-1) = snow_layer_temp(map_snow(1)+n-1) + &
          snow_layer_temp_inc(map_snow(1)+n-1)
          snow_layer_temp(map_snow(1)+n-1) = MIN(snow_layer_temp(map_snow(1)+n-1), T_freeze_h2o)
        end if

      end do

    end if

  end subroutine surf_inc_update_code

end module surf_inc_update_kernel_mod
