!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Perform IAU updates using level1 temperature increments
module level_one_temp_inc_kernel_mod

  use argument_mod,               only: arg_type,                  &
                                        GH_FIELD, GH_SCALAR,       &
                                        GH_REAL, GH_INTEGER,       &
                                        GH_READ, GH_READWRITE,     &
                                        ANY_DISCONTINUOUS_SPACE_1, &
                                        ANY_DISCONTINUOUS_SPACE_2, &
                                        ANY_DISCONTINUOUS_SPACE_3, &
                                        CELL_COLUMN
  use constants_mod,              only: i_def, r_def, l_def, r_um, rmdi
  use fs_continuity_mod,          only: WTHETA
  use jules_snow_mod,             only: dzsnow
  use kernel_mod,                 only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: level_one_temp_inc_kernel_type
    private
    type(arg_type) :: meta_args(13) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),                    & !temperature_inc
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & !soil_temperature
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & !tile_fraction
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & !snow_depth
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,   ANY_DISCONTINUOUS_SPACE_2), & !n_snow_layers
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & !tile_temperature
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_3), & !snow_layer_temp
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !n_land_tile
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !snow_lev_tile
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !sm_levels
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !nsmax
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              & !lake
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  & !T_freeze_h2o
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: level_one_temp_inc_kernel_code
  end type

  public :: level_one_temp_inc_kernel_code

contains

  !> @param[in]     nlayers           number of layers
  !> @param[in]     temperature_inc   temperature increment on theta levels
  !> @param[in,out] soil_temperature  soil temperature on soil levels
  !> @param[in]     tile_fraction     tile fraction
  !> @param[in]     snow_depth        snow depth on surface tiles
  !> @param[in]     n_snow_layers     number of snow layers on surface tiles
  !> @param[in,out] tile_temperature  tile temperature on surface tiles
  !> @param[in,out] snow_layer_temp   snow temp on snow layers and tiles
  !> @param[in]     n_land_tile       number of land tiles
  !> @param[in]     snow_lev_tile     number of snow levels and land tiles
  !> @param[in]     sm_levels         number of soil levels
  !> @param[in]     nsmax             number of snow layers in multilayer snow scheme
  !> @param[in]     lake              tile number for lake surface type
  !> @param[in]     T_freeze_h2o      freezing point of water
  !> @param[in]     ndf_wth           number of total DOFs per cell for wtheta fields
  !> @param[in]     undf_wth          number of total DOFs for wtheta fields
  !> @param[in]     map_2d            dofmap for cell for wtheta fields
  !> @param[in]     ndf_soil          number of DOFs per call for soil levels
  !> @param[in]     undf_soil         number of total DOFs for soil levels
  !> @param[in]     map_soil          Dofmap for cell for soil levels
  !> @param[in]     ndf_tile          number of DOFs per call for surface fields
  !> @param[in]     undf_tile         number of total DOFs for surface fields
  !> @param[in]     map_tile          Dofmap for cell for surface fields
  !> @param[in]     ndf_snow          number of DOFs per call for snow fields
  !> @param[in]     undf_snow         number of total DOFs for snow fields
  !> @param[in]     map_snow          Dofmap for cell for snow fields

  subroutine level_one_temp_inc_kernel_code ( nlayers,                       &
                                              temperature_inc,               &  !wth
                                              soil_temperature,              &  !soil
                                              tile_fraction,                 &  !surft
                                              snow_depth,                    &  !surft
                                              n_snow_layers,                 &  !surft
                                              tile_temperature,              &  !surft
                                              snow_layer_temp,               &  !snow
                                              n_land_tile,                   &
                                              snow_lev_tile,                 &
                                              sm_levels,                     &
                                              nsmax,                         &
                                              lake,                          &
                                              T_freeze_h2o,                  &
                                              ndf_wth, undf_wth, map_wth,    &  !wth
                                              ndf_soil, undf_soil, map_soil, &  !soil
                                              ndf_tile, undf_tile, map_tile, &  !surft
                                              ndf_snow, undf_snow, map_snow )   !snow

    implicit none

    ! arguments
    integer( kind=i_def ), intent(in) :: nlayers
    integer( kind=i_def ), intent(in) :: n_land_tile
    integer( kind=i_def ), intent(in) :: snow_lev_tile
    integer( kind=i_def ), intent(in) :: sm_levels
    integer( kind=i_def ), intent(in) :: nsmax
    integer( kind=i_def ), intent(in) :: lake

    integer( kind=i_def ), intent(in) :: ndf_wth, undf_wth
    integer( kind=i_def ), intent(in) :: map_wth(ndf_wth)
    integer( kind=i_def ), intent(in) :: ndf_soil, undf_soil
    integer( kind=i_def ), intent(in) :: map_soil(ndf_soil)
    integer( kind=i_def ), intent(in) :: ndf_tile, undf_tile
    integer( kind=i_def ), intent(in) :: map_tile(ndf_tile)
    integer( kind=i_def ), intent(in) :: ndf_snow, undf_snow
    integer( kind=i_def ), intent(in) :: map_snow(ndf_snow)

    real( kind=r_def ), intent(in)    :: T_freeze_h2o

    real( kind=r_def ), intent(in)    :: temperature_inc(undf_wth)
    real( kind=r_def ), intent(inout) :: soil_temperature(undf_soil)
    real( kind=r_def ), intent(in)    :: tile_fraction(undf_tile)
    real( kind=r_def ), intent(in)    :: snow_depth(undf_tile)
    integer( kind=i_def ), intent(in) :: n_snow_layers(undf_tile)
    real( kind=r_def ), intent(inout) :: tile_temperature(undf_tile)
    real( kind=r_def ), intent(inout) :: snow_layer_temp(undf_snow)

    ! local variables
    real(kind=r_def), parameter     :: snow_response_depth = 0.3_r_def
    real(kind=r_def), parameter     :: thin_snow_thresh = 0.04_r_def
    real(kind=r_def)                :: t_level1_inc
    real(kind=r_def)                :: total_land_frac
    real(kind=r_def)                :: land_tile_frac(n_land_tile)
    real(kind=r_def)                :: weight_soil_snow(n_land_tile)
    real(kind=r_def)                :: weight_snow_layer(nsmax)
    real(kind=r_def)                :: z_snow
    integer(kind=i_def)             :: n, i, i_snow
    logical(kind=l_def)             :: l_update_soil_temperature(n_land_tile)
    logical(kind=l_def)             :: l_update_tile_temperature(n_land_tile)
    logical(kind=l_def)             :: l_update_snow_layer_temp(snow_lev_tile)
    logical(kind=l_def)             :: no_snow(n_land_tile)
    logical(kind=l_def)             :: valid_tile(n_land_tile)

    ! Find the increment to level one temperature
    t_level1_inc = temperature_inc(map_wth(1)+1)

    ! Work out fraction of all land area for each land tile (the field
    ! tile_frac gives fraction of the cell, including sea and sea-ice tiles)
    ! Work out whether snow is present for each land tile
    ! For each land tile determine whether it is allowed to be updated - updates
    ! only allowed for non-lake tiles with > 0 tile fraction.

    ! initialisation of local variables
    total_land_frac = 0.0_r_def
    no_snow(:) = .true.
    valid_tile(:) = .false.
    l_update_tile_temperature(:) = .false.
    l_update_soil_temperature(:) = .false.
    l_update_snow_layer_temp(:) = .false.

    do n = 1, n_land_tile
      total_land_frac = total_land_frac + tile_fraction(map_tile(1)+n-1)
    end do

    if ( total_land_frac > 0.0_r_def ) then
      do n = 1, n_land_tile
        land_tile_frac(n) = tile_fraction(map_tile(1)+n-1) / total_land_frac
      end do
    else
      land_tile_frac(:) = 0.0_r_def
    end if

    do n = 1, n_land_tile

      if ( (n_snow_layers(map_tile(1)+n-1)) > 0_i_def ) no_snow(n) = .false.

      if ( (land_tile_frac(n) > 0.0_r_def) .and. (n /= lake) ) valid_tile(n) = .true.

    ! work out update switches for each field, based on the field containing valid data
      if ( (valid_tile(n)) .and. &
        (tile_temperature(map_tile(1)+n-1) > rmdi) ) then
        l_update_tile_temperature(n) = .true.
      end if

      if ( (valid_tile(n)) .and. &
        (soil_temperature(map_soil(1)) > rmdi ) ) then
        l_update_soil_temperature(n) = .true.
      end if

      do i = 1, nsmax
        i_snow = (n-1) * nsmax + i-1
        if ( (valid_tile(n)) .and. &
          (snow_layer_temp(map_snow(1)+i_snow) > rmdi) ) then
          l_update_snow_layer_temp(i_snow+1) = .true.
        end if
      end do

    end do

    ! calculate the weighting factor for level 1 temperature increments added
    ! to soil temperature under snow. The weighting ensures soil temperature is
    ! only updated if the snow is very thin. The threshold is partly historical
    ! and a proper treatment would take full account of the physical timescales
    ! in the DA interval
    do n = 1, n_land_tile

      weight_soil_snow(n) = MAX ( 0.0_r_def, 1.0_r_def - &
                            snow_depth(map_tile(1)+n-1) / thin_snow_thresh )

    end do

    ! Calculate weighting factors for each snow layer, for level 1 temperature
    ! increments added to snow layer temperature. We expect the first snow
    ! layer to be so thin that we shall always want to increment it. Typically
    ! the phase of the diurnal cycle will be reversed at a depth of about
    ! 0.3 m. The weight applied to the increment is therefore decreased
    ! smoothly from half this depth (diffusive timescale 4 times smaller) to
    ! 0 at 0.3 m.

    weight_snow_layer(1) = 1.0_r_def
    z_snow = 0.5_r_def * dzsnow(1)
    do i = 2, nsmax

      z_snow = z_snow + 0.5_r_def * ( dzsnow(i-1) + dzsnow(i) )
      weight_snow_layer(i) = ( snow_response_depth - z_snow) / &
                             ( 0.5_r_def * snow_response_depth )
      ! constrain weighting factor to be between 0 and 1
      weight_snow_layer(i) = MIN( 1.0_r_def, weight_snow_layer(i) )
      weight_snow_layer(i) = MAX( 0.0_r_def, weight_snow_layer(i) )

    end do


    ! Update each variable in turn

    do n = 1, n_land_tile

      if ( no_snow(n) ) then
     ! update tile_temperature and soil temperature (1st soil level only)
     ! soil_temperature is not a tiled field, so increment is aggregated
     ! over the non-lake land portion of the grid-box
        if ( l_update_tile_temperature(n) ) then
          tile_temperature(map_tile(1)+n-1) = tile_temperature(map_tile(1)+n-1) + &
          t_level1_inc
        end if

        if ( l_update_soil_temperature(n) ) then
          soil_temperature(map_soil(1)) = soil_temperature(map_soil(1)) + &
          land_tile_frac(n) * t_level1_inc
        end if

      else ! there is snow

        if ( (l_update_tile_temperature(n)) .and. &
          (tile_temperature(map_tile(1)+n-1) < T_freeze_h2o) ) then

          tile_temperature(map_tile(1)+n-1) = &
          MIN ( tile_temperature(map_tile(1)+n-1) + t_level1_inc, T_freeze_h2o )
        end if

        if ( (l_update_soil_temperature(n)) .and. &
          (soil_temperature(map_soil(1)) < T_freeze_h2o) ) then

          soil_temperature(map_soil(1)) = MIN ( soil_temperature(map_soil(1)) &
          + weight_soil_snow(n) * land_tile_frac(n) * t_level1_inc, T_freeze_h2o )
        end if

        do i = 1, nsmax
          i_snow = (n-1) * nsmax + i-1

          if ( l_update_snow_layer_temp(i_snow+1) ) then
            snow_layer_temp(map_snow(1)+i_snow) = &
            MIN (snow_layer_temp(map_snow(1)+i_snow) + weight_snow_layer(i) * &
                 t_level1_inc, T_freeze_h2o)
          end if

        end do

      end if ! snow or no snow

    end do

  end subroutine level_one_temp_inc_kernel_code

end module level_one_temp_inc_kernel_mod
