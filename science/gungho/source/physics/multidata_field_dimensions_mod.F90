!-------------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief The module multidata_field_dimensions_mod provides access to the
!> config variables for multidata fields (tile fields).

module multidata_field_dimensions_mod

  use constants_mod,             only: i_def, l_def
  use io_config_mod,             only: use_xios_io
  use lfric_xios_diag_mod,       only: set_axis_dimension, &
                                       get_axis_dimension, &
                                       set_zoom_axis_attr
#ifdef UM_PHYSICS
  use jules_surface_mod,              only: l_urban2t
  use jules_surface_types_config_mod, only: c3_grass,    &
                                            ice,         &
                                            urban,       &
                                            urban_canyon
#endif

  implicit none

  private

#ifdef UM_PHYSICS
      !                   1         2         3
      !          123456789012345678901234567890
      character(30), parameter :: multidata_items(33) = &
            [character(30) ::                           &
                'plant_func_types',                     &
                'sea_ice_categories',                   &
                'land_tiles',                           &
                'surface_tiles',                        &
                'surface_regrid_vars',                  &
                'snow_layers_and_tiles',                &
                'soil_levels',                          &
                'soil_levels_and_tiles',                &
                'land_tile_rad_band',                   &
                'monthly_climatology',                  &
                'boundary_layer_types',                 &
                'entrainment_levels',                   &
                'dust_divisions',                       &
                'radiation_levels',                     &
                'aero_modes',                           &
                'sw_bands_aero_modes',                  &
                'lw_bands_aero_modes',                  &
                'cloud_subcols',                        &
                'isccp_ctp_tau_bins',                   &
                'cloudsat_levels',                      &
                'csat_lvls_atb_bins',                   &
                'horizon_angles',                       &
                'horizon_aspects',                      &
                'aod_wavel',                            &
                'sw_bands',                             &
                'lw_bands',                             &
                'sw_bands_surface_tiles',               &
                'lw_bands_surface_tiles',               &
                'sw_bands_radiation_levels',            &
                'lw_bands_radiation_levels',            &
                'photolysis_pathways',                  &
                'random_seed_size',                     &
                'photol_species'                        &
      ]
#endif

  public :: get_multidata_field_dimension, sync_multidata_field_dimensions

contains

!> @brief Synchronise XIOS axis dimensions with multidata field
!> dimensions sourced from configuration
subroutine sync_multidata_field_dimensions()
      implicit none
#ifdef UM_PHYSICS
      logical(l_def), parameter :: tolerate_missing_axes = .true.
      integer(i_def) :: i
      do i=1,size(multidata_items)
            call set_axis_dimension(                                          &
                  multidata_items(i),                                         &
                  get_multidata_field_dimension(multidata_items(i)),          &
                  tolerate_missing_axes)
      end do

      ! The zoom filters in the urbanXt axis def files are something that
      ! also needs to change with the size of the land_tiles array
      ! (the begin element is n_land_tiles and n_land_tiles+1 respectively)
      call set_zoom_axis_attr(                                                &
            'surface_tiles_sea_zoom_axis',                                    &
            get_multidata_field_dimension('land_tiles'),                      &
            1,                                                                &
            tolerate_missing_axes)
      call set_zoom_axis_attr(                                                &
            'surface_tiles_sea_ice_zoom_axis',                                &
            get_multidata_field_dimension('land_tiles')+1,                    &
            1,                                                                &
            tolerate_missing_axes)
      ! Settings for some land tile zoom axes
      call set_zoom_axis_attr(                                                &
            'surface_tiles_c3_grass_zoom_axis',                               &
            (c3_grass - 1),                                                   &
            1,                                                                &
            tolerate_missing_axes)
      call set_zoom_axis_attr(                                                &
            'surface_tiles_land_ice_zoom_axis',                               &
            (ice - 1),                                                        &
            1,                                                                &
            tolerate_missing_axes)
      if ( l_urban2t ) then
        call set_zoom_axis_attr(                                                &
              'surface_tiles_urban_zoom_axis',                                  &
              (urban_canyon - 1),                                               &
              1,                                                                &
              tolerate_missing_axes)
      else
        call set_zoom_axis_attr(                                                &
              'surface_tiles_urban_zoom_axis',                                  &
              (urban - 1),                                                      &
              1,                                                                &
              tolerate_missing_axes)
      end if

#endif
end subroutine sync_multidata_field_dimensions

  !> @brief Fetches the configuration data of the number of multidata entries.
  !> @details The function get_multidata_field_dimension returns
  !> the configuration multidata size ("dimension", "tile length")
  !> of a multidata item.
  !> For instance, the multidata item of the fields "tile_fraction" and
  !> "tile_temperature" is "surface_tiles".
  !> The multidata item of the fields "ent_t_frac" and "ent_zrzi" is
  !> "entrainment_levels".
  !> @param multidata_item
  !> @return Integer containing the number of multidata entries
  function get_multidata_field_dimension(multidata_item) result (dim)

#ifdef UM_PHYSICS
    use jules_control_init_mod,  only: n_surf_tile, n_sea_ice_tile,            &
                                       soil_lev_tile, n_surf_interp, n_land_tile
    use jules_physics_init_mod,  only: snow_lev_tile
    use jules_surface_types_mod, only: npft
    use nlsizes_namelist_mod,    only: sm_levels
    use ancil_info,              only: rad_nband
    use um_physics_init_mod,     only: sw_band_mode, lw_band_mode, mode_dimen
    use socrates_init_mod,       only: n_sw_band, n_lw_band, n_pathway
    use dust_parameters_mod,     only: ndiv
    use extrusion_config_mod,    only: number_of_layers
    use cosp_config_mod,         only: n_subcol_gen
    use cosp_mod,                only: n_cloudsat_levels, n_backscatter_bins,  &
                                       n_isccp_tau_bins, n_isccp_pressure_bins
    use radiation_config_mod,    only: topography,         &
                                       topography_horizon, &
                                       n_horiz_layer,      &
                                       n_horiz_ang
    use section_choice_config_mod, &
                                 only: radiation,                              &
                                       radiation_socrates,                     &
                                       radiation_none

    use ukca_radaer_precalc,     only: npd_ukca_aod_wavel

    use aerosol_config_mod,      only: l_radaer
    use chemistry_config_mod,    only: chem_scheme, chem_scheme_strattrop
    use um_ukca_init_mod,        only: n_phot_spc
#endif

    use log_mod,                 only: log_event, LOG_LEVEL_ERROR,             &
                                       log_scratch_space

    implicit none

    character(*), intent(in) :: multidata_item

    integer(kind=i_def) :: dim

    select case (multidata_item)
#ifdef UM_PHYSICS
      case ('plant_func_types')
            dim = npft
      case ('sea_ice_categories')
            dim = n_sea_ice_tile
      case ('land_tiles')
            dim = n_land_tile
      case ('surface_tiles')
            dim = n_surf_tile
      case ('surface_regrid_vars')
            dim = n_surf_interp
      case ('snow_layers_and_tiles')
            dim = snow_lev_tile
      case ('soil_levels')
            dim = sm_levels
      case ('soil_levels_and_tiles')
            dim = soil_lev_tile
      case ('land_tile_rad_band')
            dim = n_land_tile*rad_nband
      case ('monthly_climatology')
            dim = 12
      case ('boundary_layer_types')
            dim = 7
      case ('entrainment_levels')
            dim = 3
      case ('dust_divisions')
            dim = ndiv
      case ('radiation_levels')
            dim = number_of_layers+1
      case ('aero_modes')
            dim = mode_dimen
      case ('sw_bands_aero_modes')
            dim = sw_band_mode
            if (radiation == radiation_none) dim = 1
      case ('lw_bands_aero_modes')
            dim = lw_band_mode
            if (radiation == radiation_none) dim = 1
      case ('sw_bands')
            dim = n_sw_band
            if (radiation == radiation_none) dim = 1
      case ('lw_bands')
            dim = n_lw_band
            if (radiation == radiation_none) dim = 1
      case ('lw_bands_surface_tiles')
            dim = n_lw_band*n_surf_tile
            if (radiation == radiation_none) dim = 1
      case ('sw_bands_surface_tiles')
            dim = n_sw_band*n_surf_tile
            if (radiation == radiation_none) dim = 1
      case ('sw_bands_radiation_levels')
            dim = n_sw_band*(number_of_layers+1)
            if (radiation == radiation_none) dim = 1
      case ('lw_bands_radiation_levels')
            dim = n_lw_band*(number_of_layers+1)
            if (radiation == radiation_none) dim = 1
      case ('photolysis_pathways')
            dim = n_pathway
            if (radiation == radiation_none) dim = 1
      case ('cloud_subcols')
            dim = n_subcol_gen
      case ('isccp_ctp_tau_bins')
            dim = n_isccp_tau_bins*n_isccp_pressure_bins
      case ('cloudsat_levels')
            dim = n_cloudsat_levels
      case ('csat_lvls_atb_bins')
            dim = n_cloudsat_levels*n_backscatter_bins
      case ('horizon_angles')
            if (radiation == radiation_socrates .and. &
                  topography == topography_horizon) then
                  dim = n_horiz_layer*n_horiz_ang
            else
                  dim = 1
            end if
      case ('horizon_aspects')
            if (radiation == radiation_socrates .and. &
                  topography == topography_horizon) then
                  dim = n_horiz_ang
            else
                  dim = 1
            end if
      case('aod_wavel')

           if ( l_radaer ) then
                  dim = npd_ukca_aod_wavel
           else
                  dim = 1
           end if

      case('random_seed_size')
           call random_seed(size=dim)
      case ('photol_species')
            if (chem_scheme == chem_scheme_strattrop) then
               dim = n_phot_spc
            else
               dim = 1
            end if
      case ('')
            dim = 1 ! ordinary (non-multidata) field
#endif
      case default
            if (use_xios_io) then
                  ! attempt to get it from XIOS metadata
                  dim = get_axis_dimension(multidata_item)
            else
                  dim = 1 ! silence compiler warning
                  write(log_scratch_space, '(A, A)')                          &
                        'Unexpected multidata item: ', multidata_item
                        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
            end if
    end select

  end function get_multidata_field_dimension

end module multidata_field_dimensions_mod
