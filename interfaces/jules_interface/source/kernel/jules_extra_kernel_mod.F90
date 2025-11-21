!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to JULES surf_couple_extra.
!>
module jules_extra_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_WRITE,         &
                                      GH_READWRITE, GH_INTEGER,  &
                                      ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2, &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      ANY_DISCONTINUOUS_SPACE_4, &
                                      ANY_DISCONTINUOUS_SPACE_5, &
                                      DOMAIN
  use constants_mod,           only : i_def, i_um, r_def, r_um
  use kernel_mod,              only : kernel_type
  use empty_data_mod,          only : empty_real_data

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: jules_extra_kernel_type
    private
    type(arg_type) :: meta_args(58) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! ls_rain
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! conv_rain
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! ls_snow
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! conv_snow
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! ls_graup
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! lsca_2d
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! cca_2d
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! leaf_area_index
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! canopy_height
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! snow_unload_rate
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_wilt
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_crit
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_sat
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_cond_sat
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_thermal_cap
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_thermal_cond
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_suction_sat
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! clapp_horn_b
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_roughness
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! mean_topog_index
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! a_sat_frac
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! c_sat_frac
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! a_wet_frac
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! c_wet_frac
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! net_prim_prod
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! snowice_sublimation
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! surf_heat_flux
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! canopy_evap
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_4), & ! water_extraction
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! thermal_cond_wet_soil
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! urbztm
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! soil_temperature
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! soil_moisture
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! unfrozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! frozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! canopy_water
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_snow_mass
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_snow_rgrain
         arg_type(GH_FIELD, GH_INTEGER,GH_READWRITE,ANY_DISCONTINUOUS_SPACE_2),& ! n_snow_layers
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! snow_depth
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! snow_under_canopy
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! snowpack_density
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_thickness
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_ice_mass
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_liq_mass
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_temp
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_rgrain
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! snowice_melt
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! soil_sat_frac
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! water_table
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! wetness_under_soil
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! surface_runoff
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sub_surface_runoff
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! soil_moisture_content
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! grid_snow_mass
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2)  & ! throughfall
        /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: jules_extra_code
  end type

  public :: jules_extra_code

contains

  !> @brief Interface to JULES surf_couple_extra
  !> @details JULES surf_couple_extra calculates the surface and soil fluxes
  !> and stores of Water (via the hydrol{ogy}, snow and river_control
  !> modules) and Carbon (via triffid and inferno/fire modules). Note: only the
  !> hydrol and snow components are implemented so far.

  !> @param[in]     nlayers                Number of layers
  !> @param[in]     ls_rain                Large-scale rainfall rate (kg m-2 s-1)
  !> @param[in]     conv_rain              Convective rainfall rate (kg m-2 s-1)
  !> @param[in]     ls_snow                Large-scale snowfall rate (kg m-2 s-1)
  !> @param[in]     conv_snow              Convective snowfall rate (kg m-2 s-1)
  !> @param[in]     ls_graup               Large-scale graupelfall rate (kg m-2 s-1)
  !> @param[in]     lsca_2d                Large-scale cloud amout (2d)
  !> @param[in]     cca_2d                 Convective cloud amout (2d) with no anvil
  !> @param[in]     tile_fraction          Surface tile fractions
  !> @param[in]     leaf_area_index        Leaf Area Index
  !> @param[in]     canopy_height          Canopy height (m)
  !> @param[in]     snow_unload_rate       Unloading of snow from PFTs by wind
  !> @param[in]     soil_moist_wilt        Volumetric soil moist at wilting pt
  !> @param[in]     soil_moist_crit        Volumetric soil moist at critical pt
  !> @param[in]     soil_moist_sat         Volumetric soil moist at saturation
  !> @param[in]     soil_cond_sat          Saturated soil thermal conductivity (kg m-2 s-1)
  !> @param[in]     soil_thermal_cap       Soil thermal capacity (J m-3 K-1)
  !> @param[in]     soil_thermal_cond      Soil thermal conductivity (W m-1 K-1)
  !> @param[in]     soil_suction_sat       Saturated soil water suction (m)
  !> @param[in]     clapp_horn_b           Clapp and Hornberger b coefficient
  !> @param[in]     soil_roughness         Bare soil surface roughness length (m)
  !> @param[in]     mean_topog_index       Mean topographic index
  !> @param[in]     a_sat_frac             a gridbox saturated fraction
  !> @param[in]     c_sat_frac             c gridbox saturated fraction
  !> @param[in]     a_wet_frac             a gridbox wet fraction
  !> @param[in]     c_wet_frac             c gridbox wet fraction
  !> @param[in]     tile_temperature       Surface tile temperatures (K)
  !> @param[in]     net_prim_prod          Net Primary Productivity (kg m-2 s-1)
  !> @param[in]     snowice_sublimation    Sublimation of snow and ice (kg m-2 s-1)
  !> @param[in]     surf_heat_flux         Surface heat flux (W m-2)
  !> @param[in]     canopy_evap            Canopy evaporation from land tiles (kg m-2 s-1)
  !> @param[in]     water_extraction       Extraction of water from each soil layer (kg m-2 s-1)
  !> @param[in]     thermal_cond_wet_soil  Thermal conductivity of soil (W m-1 K-1)
  !> @param[in]     urbztm                 Urban effective roughness length
  !> @param[in,out] soil_temperature       Soil temperature (K)
  !> @param[in,out] soil_moisture          Soil moisture content (kg m-2)
  !> @param[in,out] unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in,out] frozen_soil_moisture   Frozen soil moisture proportion
  !> @param[in,out] canopy_water           Canopy water on each tile (kg m-2)
  !> @param[in,out] tile_snow_mass         Snow mass on tiles (kg m-2)
  !> @param[in,out] tile_snow_rgrain       Snow grain radius on tiles (microns)
  !> @param[in,out] n_snow_layers          Number of snow layers on tiles
  !> @param[in,out] snow_depth             Snow depth on tiles (m)
  !> @param[in,out] snow_under_canopy      Amount of snow under canopy (kg m-2)
  !> @param[in,out] snowpack_density       Density of snow on ground (kg m-3)
  !> @param[in,out] snow_layer_thickness   Thickness of snow layers (m)
  !> @param[in,out] snow_layer_ice_mass    Mass of ice in snow layers (kg m-2)
  !> @param[in,out] snow_layer_liq_mass    Mass of liquid in snow layers (kg m-2)
  !> @param[in,out] snow_layer_temp        Temperature of snow layer (K)
  !> @param[in,out] snow_layer_rgrain      Grain radius of snow layer (microns)
  !> @param[in,out] snowice_melt           Surface, canopy and sea ice, snow and ice melt rate (kg m-2 s-1)
  !> @param[in,out] soil_sat_frac          Soil saturated fraction
  !> @param[in,out] water_table            Water table depth (m)
  !> @param[in,out] wetness_under_soil     Soil wetness below soil column
  !> @param[in,out] surface_runoff         Runoff from surface
  !> @param[in,out] sub_surface_runoff     Runoff from sub-surface
  !> @param[in,out] soil_moisture_content  Soil moisture content of soil column
  !> @param[in,out] grid_snow_mass         Gridbox total snow mass (canopy + under canopy)
  !> @param[in,out] throughfall            Throughfall from land tiles
  !> @param[in]     ndf_2d                 Total DOFs per cell for 2D fields
  !> @param[in]     undf_2d                Unique DOFs per cell for 2D fields
  !> @param[in]     map_2d                 DOFmap for cells for 2D fields
  !> @param[in]     ndf_tile               Total DOFs per cell for surface tiles
  !> @param[in]     undf_tile              Unique DOFs per cell for surface tiles
  !> @param[in]     map_tile               DOFmap for cells for surface tiles
  !> @param[in]     ndf_pft                Total DOFs per cell for plant types
  !> @param[in]     undf_pft               Unique DOFs per cell for plant types
  !> @param[in]     map_pft                DOFmap for cells for plant types
  !> @param[in]     ndf_soil               Total DOFs per cell for soil levels
  !> @param[in]     undf_soil              Unique DOFs per cell for soil levels
  !> @param[in]     map_soil               DOFmap for cells for soil levels
  !> @param[in]     ndf_snow               Total DOFs per cell for snow layers
  !> @param[in]     undf_snow              Unique DOFs per cell for snow layers
  !> @param[in]     map_snow               DOFmap for cells for snow layers
  subroutine jules_extra_code(             &
               nlayers, seg_len,           &
               ls_rain,                    &
               conv_rain,                  &
               ls_snow,                    &
               conv_snow,                  &
               ls_graup,                   &
               lsca_2d,                    &
               cca_2d,                     &
               tile_fraction,              &
               leaf_area_index,            &
               canopy_height,              &
               snow_unload_rate,           &
               soil_moist_wilt,            &
               soil_moist_crit,            &
               soil_moist_sat,             &
               soil_cond_sat,              &
               soil_thermal_cap,           &
               soil_thermal_cond,          &
               soil_suction_sat,           &
               clapp_horn_b,               &
               soil_roughness,             &
               mean_topog_index,           &
               a_sat_frac,                 &
               c_sat_frac,                 &
               a_wet_frac,                 &
               c_wet_frac,                 &
               tile_temperature,           &
               net_prim_prod,              &
               snowice_sublimation,        &
               surf_heat_flux,             &
               canopy_evap,                &
               water_extraction,           &
               thermal_cond_wet_soil,      &
               urbztm,                     &
               soil_temperature,           &
               soil_moisture,              &
               unfrozen_soil_moisture,     &
               frozen_soil_moisture,       &
               canopy_water,               &
               tile_snow_mass,             &
               tile_snow_rgrain,           &
               n_snow_layers,              &
               snow_depth,                 &
               snow_under_canopy,          &
               snowpack_density,           &
               snow_layer_thickness,       &
               snow_layer_ice_mass,        &
               snow_layer_liq_mass,        &
               snow_layer_temp,            &
               snow_layer_rgrain,          &
               snowice_melt,               &
               soil_sat_frac,              &
               water_table,                &
               wetness_under_soil,         &
               surface_runoff,             &
               sub_surface_runoff,         &
               soil_moisture_content,      &
               grid_snow_mass,             &
               throughfall,                &
               ndf_2d,                     &
               undf_2d,                    &
               map_2d,                     &
               ndf_tile,                   &
               undf_tile,                  &
               map_tile,                   &
               ndf_pft,                    &
               undf_pft,                   &
               map_pft,                    &
               ndf_soil,                   &
               undf_soil,                  &
               map_soil,                   &
               ndf_snow,                   &
               undf_snow,                  &
               map_snow                    &
                            )

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use jules_control_init_mod,     only: nsurft => n_land_tile
    use jules_physics_init_mod,     only: decrease_sath_cond

    ! Module imports for surf_couple_extra JULESvn5.4
    use ancil_info,               only: nsoilt, dim_cslayer, rad_nband,        &
                                        dim_soil_n_pool, nmasst
    use atm_step_local,           only: dim_cs1
    use cderived_mod,             only: delta_lambda, delta_phi
    use dust_parameters_mod,      only: ndiv
    use jules_hydrology_mod,      only: l_hydrology
    use jules_surface_types_mod,  only: npft, ntype
    use jules_snow_mod,           only: nsmax
    use jules_soil_mod,           only: ns_deep, l_bedrock
    use jules_soil_biogeochem_mod, only: dim_ch4layer, soil_bgc_model,         &
                                         soil_model_ecosse, l_layeredc
    use jules_snow_mod,           only: nsmax, cansnowtile
    use jules_deposition_mod,     only: l_deposition
    use jules_sea_seaice_mod,     only: nice, nice_use
    use jules_deposition_mod,     only: l_deposition
    use jules_surface_mod,        only: l_urban2t, l_flake_model
    use jules_urban_mod,          only: l_moruses
    use jules_vegetation_mod,     only: l_crop, l_triffid,                     &
                                        l_phenol, l_use_pft_psi, can_rad_mod,  &
                                        l_acclim, l_sugar, l_red
    use theta_field_sizes,        only: t_i_length, t_j_length, u_i_length,    &
                                        u_j_length, v_i_length, v_j_length
    use jules_surface_types_mod,  only: ncpft, nnpft
    use jules_irrig_mod,          only: irr_crop, irr_crop_doell
    use atm_fields_bounds_mod,    only: pdims_s, pdims
    use jules_radiation_mod,      only: l_albedo_obs
    use jules_water_resources_mod,                                             &
                                  only: l_have_groundwater,                    &
                                        l_have_surface_water, l_water_domestic,&
                                        l_water_industry, l_water_irrigation,  &
                                        l_water_livestock, l_water_resources,  &
                                        l_water_transfers, n_sw_source,        &
                                        nwater_use
    use jules_water_tracers_mod,  only: l_wtrac_jls, n_wtrac_jls, n_evap_srce

    use crop_vars_mod,            only: crop_vars_type, crop_vars_data_type,   &
                                        crop_vars_alloc, crop_vars_assoc, &
                                        crop_vars_dealloc, crop_vars_nullify
    use prognostics,              only: progs_data_type, progs_type,           &
                                        prognostics_alloc, prognostics_assoc,  &
                                        prognostics_nullify, prognostics_dealloc
    use fire_vars_mod,            only: fire_vars_type, fire_vars_data_type,   &
                                        fire_vars_alloc, fire_vars_assoc, &
                                        fire_vars_dealloc, fire_vars_nullify
    use jules_vars_mod,           only: jules_vars_type, jules_vars_data_type, &
                                        jules_vars_alloc, jules_vars_assoc,    &
                                        jules_vars_dealloc, jules_vars_nullify
    use p_s_parms,                only: psparms_type, psparms_data_type,       &
                                        psparms_alloc, psparms_assoc, &
                                        psparms_dealloc, psparms_nullify
    use top_pdm,                  only: top_pdm_type, top_pdm_data_type,       &
                                        top_pdm_assoc, top_pdm_alloc, &
                                        top_pdm_dealloc, top_pdm_nullify
    use trif_vars_mod,            only: trif_vars_type, trif_vars_data_type,   &
                                        trif_vars_alloc, trif_vars_assoc, &
                                        trif_vars_nullify, trif_vars_dealloc
    use soil_ecosse_vars_mod,     only: soil_ecosse_vars_type,                 &
                                        soil_ecosse_vars_data_type,            &
                                        soil_ecosse_vars_alloc,                &
                                        soil_ecosse_vars_assoc, &
                                        soil_ecosse_vars_dealloc, &
                                        soil_ecosse_vars_nullify
    use urban_param_mod,          only: urban_param_data_type,                 &
                                        urban_param_type, &
                                        urban_param_alloc, urban_param_assoc,  &
                                        urban_param_dealloc, urban_param_nullify
    use trifctl,                  only: trifctl_type, trifctl_data_type,       &
                                        trifctl_alloc, trifctl_assoc,          &
                                        trifctl_nullify, trifctl_dealloc
    use lake_mod,                 only: lake_type, lake_data_type,             &
                                        lake_alloc, lake_assoc,                &
                                        lake_nullify, lake_dealloc
    use ancil_info,               only: ainfo_type, ainfo_data_type,           &
                                        ancil_info_assoc, ancil_info_alloc,    &
                                        ancil_info_dealloc, ancil_info_nullify
    use jules_forcing_mod,        only: forcing_type, forcing_data_type,       &
                                        forcing_assoc, forcing_alloc,          &
                                        forcing_nullify, forcing_dealloc
    use fluxes_mod,               only: fluxes_type, fluxes_data_type,         &
                                        fluxes_alloc, fluxes_assoc,            &
                                        fluxes_nullify, fluxes_dealloc
    use jules_rivers_mod,         only: rivers_type, rivers_data_type,         &
                                        rivers_assoc, jules_rivers_alloc,      &
                                        rivers_nullify, rivers_dealloc
    use cable_fields_mod,         only: work_vars_cbl
    use jules_chemvars_mod,       only: chemvars_type, chemvars_data_type,     &
                                        chemvars_alloc, chemvars_assoc,        &
                                        chemvars_nullify, chemvars_dealloc
    use water_resources_vars_mod, only: water_resources_type,                  &
                                        water_resources_data_type,             &
                                        water_resources_alloc,                 &
                                        water_resources_assoc,                 &
                                        water_resources_nullify,               &
                                        water_resources_dealloc
    use jules_wtrac_type_mod,     only: jls_wtrac_type, jls_wtrac_data_type,   &
                                        wtrac_jls_assoc, wtrac_jls_alloc,      &
                                        wtrac_jls_nullify, wtrac_jls_dealloc
    use coastal,                  only: coastal_type

    use nlsizes_namelist_mod, only: sm_levels, ntiles, bl_levels
    use UM_ParCore, only: nproc

    ! Jules related subroutines
    use sparm_mod, only             : sparm
    use tilepts_mod, only           : tilepts
    use surf_couple_extra_mod, only : surf_couple_extra
    use infiltration_rate_mod, only : infiltration_rate

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d, ndf_tile, undf_tile,  &
                                       ndf_pft, undf_pft, ndf_soil,           &
                                       undf_soil, ndf_snow, undf_snow

    integer(kind=i_def), dimension(ndf_2d, seg_len),    intent(in) :: map_2d
    integer(kind=i_def), dimension(ndf_tile, seg_len),  intent(in) :: map_tile
    integer(kind=i_def), dimension(ndf_pft, seg_len),   intent(in) :: map_pft
    integer(kind=i_def), dimension(ndf_soil, seg_len),  intent(in) :: map_soil
    integer(kind=i_def), dimension(ndf_snow, seg_len),  intent(in) :: map_snow

    real(kind=r_def), dimension(undf_2d), intent(in) :: ls_rain, conv_rain,   &
                                                        ls_snow, conv_snow,   &
                                                        ls_graup, lsca_2d,    &
                                                        cca_2d

    real(kind=r_def), intent(in)    :: tile_fraction(undf_tile)
    real(kind=r_def), intent(in)    :: snowice_sublimation(undf_tile)
    real(kind=r_def), intent(in)    :: surf_heat_flux(undf_tile)
    real(kind=r_def), intent(in)    :: canopy_evap(undf_tile)
    real(kind=r_def), intent(in)    :: tile_temperature(undf_tile)

    real(kind=r_def), intent(in)    :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in)    :: canopy_height(undf_pft)
    real(kind=r_def), intent(in)    :: snow_unload_rate(undf_pft)

    real(kind=r_def), intent(in)    :: soil_moist_wilt(undf_2d)
    real(kind=r_def), intent(in)    :: soil_moist_crit(undf_2d)
    real(kind=r_def), intent(in)    :: soil_moist_sat(undf_2d)
    real(kind=r_def), intent(in)    :: soil_cond_sat(undf_2d)
    real(kind=r_def), intent(in)    :: soil_thermal_cap(undf_2d)
    real(kind=r_def), intent(in)    :: soil_suction_sat(undf_2d)
    real(kind=r_def), intent(in)    :: clapp_horn_b(undf_2d)
    real(kind=r_def), intent(in)    :: water_extraction(undf_soil)

    real(kind=r_def), intent(in)    :: soil_thermal_cond(undf_2d)
    real(kind=r_def), intent(in)    :: soil_roughness(undf_2d)
    real(kind=r_def), intent(in)    :: mean_topog_index(undf_2d)
    real(kind=r_def), intent(in)    :: a_sat_frac(undf_2d)
    real(kind=r_def), intent(in)    :: c_sat_frac(undf_2d)
    real(kind=r_def), intent(in)    :: a_wet_frac(undf_2d)
    real(kind=r_def), intent(in)    :: c_wet_frac(undf_2d)
    real(kind=r_def), intent(in)    :: net_prim_prod(undf_2d)
    real(kind=r_def), intent(in)    :: thermal_cond_wet_soil(undf_2d)
    real(kind=r_def), intent(in)    :: urbztm(undf_2d)

    real(kind=r_def), intent(inout) :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_snow_mass(undf_tile)
    real(kind=r_def), intent(inout) :: tile_snow_rgrain(undf_tile)
    integer(kind=i_def), intent(inout) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(inout) :: snow_depth(undf_tile)
    real(kind=r_def), intent(inout) :: snow_under_canopy(undf_tile)
    real(kind=r_def), intent(inout) :: snowpack_density(undf_tile)
    real(kind=r_def), intent(inout) :: snowice_melt(undf_tile)

    real(kind=r_def), intent(inout) :: snow_layer_thickness(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_ice_mass(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_liq_mass(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_temp(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_rgrain(undf_snow)

    real(kind=r_def), intent(inout) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(inout) :: soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: frozen_soil_moisture(undf_soil)

    real(kind=r_def), intent(inout) :: soil_sat_frac(undf_2d)
    real(kind=r_def), intent(inout) :: water_table(undf_2d)
    real(kind=r_def), intent(inout) :: wetness_under_soil(undf_2d)
    real(kind=r_def), intent(inout) :: surface_runoff(undf_2d)
    real(kind=r_def), intent(inout) :: sub_surface_runoff(undf_2d)

    real(kind=r_def), pointer, intent(inout) :: soil_moisture_content(:)
    real(kind=r_def), pointer, intent(inout) :: grid_snow_mass(:)
    real(kind=r_def), pointer, intent(inout) :: throughfall(:)

    ! Local variables for the kernel
    integer(kind=i_def) :: i, j, n, i_snow, m, l

!------------------------------------------------------------------------------
    ! JULES surf_couple_extra subroutine arguments declared using JULESvn5.4
    !       order and naming

    ! Integer parameters
    integer(i_um), parameter :: river_row_length = 1, river_rows = 1,         &
                                aocpl_row_length = 1, aocpl_p_rows = 1

    ! Integer indices (module intent=in)
    integer(i_um) :: a_step, g_p_field, g_r_field, global_row_length,         &
         global_rows, global_river_row_length, global_river_rows,             &
         land_pts, lice_pts, soil_pts

    ! Logical (module intent=in)
    logical :: smlt, stf_sub_surf_roff

    logical, dimension(seg_len, 1) :: land_sea_mask

    ! Real variables (module intent=in)
    ! Driving data
    real(r_um), dimension(seg_len, 1) :: ls_graup_ij,                         &
         u_1_ij, v_1_ij, cca_2d_ij, soil_clay_ij, cos_theta_latitude,         &
         flash_rate_ancil, pop_den_ancil, wealth_index_ancil, flandg, rho_star

    ! State
    real(r_um), dimension(:,:), allocatable :: u_s_std_surft

    real(r_um), dimension(:,:), allocatable :: fexp_soilt, gamtot_soilt,      &
         ti_mean_soilt, ti_sig_soilt, a_fsat_soilt, c_fsat_soilt,             &
         a_fwet_soilt, c_fwet_soilt

    real(r_um), dimension(:), allocatable :: npp_gb, frac_agr_gb

    ! River routing
    real(r_um), dimension(aocpl_row_length) :: xpa
    real(r_um), dimension(0:aocpl_row_length) :: xua
    real(r_um), dimension(aocpl_row_length+1) :: xva
    real(r_um), dimension(aocpl_p_rows) :: ypa, yua
    real(r_um), dimension(0:aocpl_p_rows) :: yva
    real(r_um), dimension(river_row_length, river_rows) :: trivdir, trivseq
    real(r_um), dimension(seg_len, 1) :: r_area, slope, flowobs1,       &
         r_inext, r_jnext, r_land


    ! Integers (module intent = in out)
    integer(i_um) :: a_steps_since_riv, asteps_since_triffid

    ! Real variables (module intent = in out)
    real(r_um), dimension(seg_len, 1) :: substore, surfstore, flowin,   &
         bflowin, acc_lake_evap

    real(r_um), dimension(river_row_length, river_rows) :: twatstor


    real(r_um), dimension(:), allocatable :: dhf_surf_minus_soil,             &
         ls_rainfrac_gb, tot_surf_runoff, tot_sub_runoff, inlandout_atm_gb

    real(r_um), dimension(:,:), allocatable :: hcons_soilt, fsat_soilt,       &
         fwetl_soilt, zw_soilt, sthzw_soilt

    !-----------------------------------------------------------------------
    ! JULES Types
    !-----------------------------------------------------------------------
    type(crop_vars_type) :: crop_vars
    type(crop_vars_data_type) :: crop_vars_data
    type(progs_type) :: progs
    type(progs_data_type) :: progs_data
    type(fire_vars_type) :: fire_vars
    type(fire_vars_data_type) :: fire_vars_data
    type(jules_vars_type) :: jules_vars
    type(jules_vars_data_type) :: jules_vars_data
    type(psparms_type) :: psparms
    type(psparms_data_type) :: psparms_data
    type(top_pdm_type) :: toppdm
    type(top_pdm_data_type) :: top_pdm_data
    type(trif_vars_type) :: trif_vars
    type(trif_vars_data_type) :: trif_vars_data
    type(soil_ecosse_vars_type) :: soilecosse
    type(soil_ecosse_vars_data_type) :: soil_ecosse_vars_data
    type(urban_param_type) :: urban_param
    type(urban_param_data_type) :: urban_param_data
    type(trifctl_type) :: trifctltype
    type(trifctl_data_type) :: trifctl_data
    type(lake_type) :: lake_vars
    type(lake_data_type) :: lake_data
    type(ainfo_type) :: ainfo
    type(ainfo_data_type) :: ainfo_data
    type(forcing_type) :: forcing
    type(forcing_data_type) :: forcing_data
    type(fluxes_type) :: fluxes
    type(fluxes_data_type) :: fluxes_data
    type(rivers_type) :: rivers
    type(rivers_data_type) :: rivers_data

    ! Variables for dry deposition
    type(chemvars_type) :: chemvars
    type(chemvars_data_type) :: chemvars_data
    type(water_resources_type) :: water_resources
    type(water_resources_data_type) :: water_resources_data
    type(jls_wtrac_type)       :: wtrac_jls
    type(jls_wtrac_data_type)  :: wtrac_jls_data
    type( coastal_type ) :: coast
    integer(i_um) :: ndry_dep_species  ! Dummy variable for now

    !-----------------------------------------------------------------------
    ! Initialisation of JULES data and pointer types
    !-----------------------------------------------------------------------
    ! Land tile fractions
    land_pts = 0
    do i = 1, seg_len
      flandg(i,1) = 0.0_r_um
      do n = 1, nsurft
        flandg(i,1) = flandg(i,1) + real(tile_fraction(map_tile(1,i)+n-1), r_um)
      end do
      flandg(i,1) = min(flandg(i,1), 1.0_r_um)
    end do

    do i = 1, seg_len
      if (flandg(i,1) > 0.0_r_um) then
        land_pts = land_pts + 1
        land_sea_mask(i,1) = .true.
      else
        land_sea_mask(i,1) = .false.
      end if
    end do

    if (land_pts == 0) then
      ! If there's no land, we can just exit here
      return
    end if

    call crop_vars_alloc(land_pts, t_i_length, t_j_length,                    &
                     nsurft, ncpft,nsoilt, sm_levels, l_crop, irr_crop,       &
                     irr_crop_doell, crop_vars_data)

    call crop_vars_assoc(crop_vars, crop_vars_data)

    call prognostics_alloc(land_pts, t_i_length, t_j_length,                  &
                      nsurft, npft, nsoilt, sm_levels, ns_deep, nsmax,        &
                      dim_cslayer, dim_cs1, dim_ch4layer,                     &
                      nice, nice_use, soil_bgc_model, soil_model_ecosse,      &
                      l_layeredc, l_triffid, l_phenol, l_bedrock, l_red,      &
                      nmasst, nnpft, l_acclim, l_sugar, progs_data)
    call prognostics_assoc(progs,progs_data)

    call fire_vars_alloc(land_pts,npft, fire_vars_data)
    call fire_vars_assoc(fire_vars, fire_vars_data)

    call jules_vars_alloc(land_pts,ntype,nsurft,rad_nband,nsoilt,sm_levels,   &
                t_i_length, t_j_length, npft, bl_levels, pdims_s, pdims,      &
                l_albedo_obs, cansnowtile, l_deposition,                      &
                jules_vars_data)
    call jules_vars_assoc(jules_vars,jules_vars_data)

    if (can_rad_mod == 6) then
      jules_vars%diff_frac = 0.4_r_um
    end if

    call psparms_alloc(land_pts,t_i_length,t_j_length,                        &
                     nsoilt,sm_levels,dim_cslayer,nsurft,npft,                &
                     soil_bgc_model,soil_model_ecosse,l_use_pft_psi,          &
                     psparms_data)
    call psparms_assoc(psparms, psparms_data)

    call top_pdm_alloc(land_pts,nsoilt, top_pdm_data)
    call top_pdm_assoc(toppdm, top_pdm_data)

    call trif_vars_alloc(land_pts,                                            &
                     npft,dim_cslayer,nsoilt,dim_cs1,                         &
                     l_triffid, l_phenol, trif_vars_data)
    call trif_vars_assoc(trif_vars, trif_vars_data)

    call soil_ecosse_vars_alloc(land_pts,                                     &
                            nsoilt,dim_cslayer,dim_soil_n_pool,sm_levels,     &
                            soil_bgc_model,soil_model_ecosse,                 &
                            soil_ecosse_vars_data)
    call soil_ecosse_vars_assoc(soilecosse, soil_ecosse_vars_data)

    call urban_param_alloc(land_pts, l_urban2t, l_moruses, urban_param_data)
    call urban_param_assoc(urban_param, urban_param_data)

    call trifctl_alloc(land_pts,                                              &
                   npft,dim_cslayer,dim_cs1,nsoilt,trifctl_data)
    call trifctl_assoc(trifctltype, trifctl_data)

    call lake_alloc(land_pts, l_flake_model, lake_data)
    call lake_assoc(lake_vars, lake_data)

    call ancil_info_alloc(land_pts,t_i_length,t_j_length,                     &
                      nice,nsoilt,ntype,                                      &
                      ainfo_data)
    call ancil_info_assoc(ainfo, ainfo_data)

    call forcing_alloc(t_i_length,t_j_length,u_i_length, u_j_length,          &
                       v_i_length, v_j_length, forcing_data)
    call forcing_assoc(forcing, forcing_data)

    call fluxes_alloc(land_pts, t_i_length, t_j_length,                       &
                      nsurft, npft, nsoilt, sm_levels,                        &
                      nice, nice_use,                                         &
                      fluxes_data)
    call fluxes_assoc(fluxes, fluxes_data)

    call jules_rivers_alloc(land_pts, t_i_length, t_j_length, rivers_data)
    call rivers_assoc(rivers,rivers_data)

    ! Chemvars for Dry deposition
    ndry_dep_species = 1
    call chemvars_alloc(land_pts, t_i_length, t_j_length, npft, ntype,        &
                        l_deposition, ndry_dep_species, chemvars_data)
    call chemvars_assoc(chemvars, chemvars_data)

    call water_resources_alloc(land_pts, n_sw_source, nwater_use,             &
                               l_have_groundwater, l_have_surface_water,      &
                               l_water_domestic, l_water_industry,            &
                               l_water_irrigation, l_water_livestock,         &
                               l_water_resources, l_water_transfers,          &
                               water_resources_data)
    call water_resources_assoc(water_resources,water_resources_data)

    call wtrac_jls_alloc(land_pts, t_i_length, t_j_length, nsurft, nsoilt,     &
                         sm_levels, nsmax, nice_use, n_wtrac_jls, n_evap_srce, &
                         river_row_length, river_rows, l_wtrac_jls,            &
                         wtrac_jls_data)
    call wtrac_jls_assoc(wtrac_jls, wtrac_jls_data)

    !-------------------------------------------------------------------
    l = 0
    do i = 1, seg_len
      if (flandg(i,1) > 0.0_r_um) then
        l = l+1
        ainfo%land_index(l) = i
      end if
    end do

    ! Data from 2D fields
    do i = 1, seg_len
      forcing%ls_rain_ij(i,1) = ls_rain(map_2d(1,i))   ! Large-scale rainfall rate
      forcing%con_rain_ij(i,1) = conv_rain(map_2d(1,i)) ! Convective rainfall rate
      forcing%ls_snow_ij(i,1) = ls_snow(map_2d(1,i))   ! Large-scale snowfall rate
      forcing%con_snow_ij(i,1) = conv_snow(map_2d(1,i)) ! Convective snowfallfall rate
      cca_2d_ij(i,1)   = cca_2d(map_2d(1,i))    ! Convective cloud amount
      ls_graup_ij(i,1) = ls_graup(map_2d(1,i)) ! Large-scale graupelfall rate
    end do

    allocate(ls_rainfrac_gb(land_pts))
    do l = 1, land_pts
      ls_rainfrac_gb(l) = lsca_2d(map_2d(1,ainfo%land_index(l)))  ! Large-scale cloud amount
    end do

    !----------------------------------------------------------------------
    ! Surface fields as needed by Jules
    !----------------------------------------------------------------------

    ! Ancillaries:
    ! Land tile fractions (frac_surft)
    do l = 1, land_pts
      do n = 1, nsurft
        ainfo%frac_surft(l, n) = real(tile_fraction(map_tile(1,ainfo%land_index(l))+n-1), r_um) &
             / flandg(ainfo%land_index(l), 1)
        fluxes%ei_surft(l, n) = real(snowice_sublimation(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        fluxes%surf_htf_surft(l, n) = real(surf_heat_flux(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        fluxes%ecan_surft(l, n) = real(canopy_evap(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      end do
    end do

    ! Logical flags controlling diagnostic calculations
    smlt = .false.
    stf_sub_surf_roff = .true.

    ! Set type_pts and type_index
    call tilepts(land_pts, ainfo%frac_surft, ainfo%surft_pts,                  &
         ainfo%surft_index, ainfo%l_lice_point, ainfo%l_lice_surft)

    ! Vegetation prognostics (for TRIFFID)
    do l = 1, land_pts
      do n = 1, npft
        ! Leaf area index
        progs%lai_pft(l, n) = real(leaf_area_index(map_pft(1,ainfo%land_index(l))+n-1), r_um)
        ! Canopy height
        progs%canht_pft(l, n) = real(canopy_height(map_pft(1,ainfo%land_index(l))+n-1), r_um)
        ! Unloading rate of snow from plant functional types
        jules_vars%unload_backgrnd_pft(l, n) = real(snow_unload_rate(map_pft(1,ainfo%land_index(l))+n-1), r_um)
      end do
    end do

    ! Soil roughness
    do l = 1, land_pts
      psparms%z0m_soil_gb(l) = real(soil_roughness(map_2d(1,ainfo%land_index(l))), r_um)
    end do
    if (l_urban2t) then
      do l = 1, land_pts
        urban_param%ztm_gb(l)  = real(urbztm(map_2d(1,ainfo%land_index(l))), r_um)
      end do
    end if

    ! Get catch_snow_surft and catch_surft from call to sparm
    call sparm(land_pts, nsurft, ainfo%surft_pts, ainfo%surft_index,           &
               ainfo%frac_surft, progs%canht_pft, progs%lai_pft,               &
               psparms%z0m_soil_gb, psparms%catch_snow_surft,                  &
               psparms%catch_surft, psparms%z0_surft, psparms%z0h_bare_surft,  &
               urban_param%ztm_gb)

    ! Decrease in saturated conductivity of soil with depth
    allocate(fexp_soilt(land_pts, nsoilt))
    do n = 1, nsoilt
      do l = 1, land_pts
        fexp_soilt(l,n) = decrease_sath_cond
      end do
    end do

    allocate(ti_mean_soilt(land_pts, nsoilt))
    allocate(a_fsat_soilt(land_pts, nsoilt))
    allocate(c_fsat_soilt(land_pts, nsoilt))
    allocate(a_fwet_soilt(land_pts, nsoilt))
    allocate(c_fwet_soilt(land_pts, nsoilt))
    allocate(npp_gb(land_pts))
    do l = 1, land_pts
      ! Mean grid box topographic index
      ti_mean_soilt(l,1) = real(mean_topog_index(map_2d(1,ainfo%land_index(l))), r_um)
      ! a saturated soil fraction
      a_fsat_soilt(l,1) = real(a_sat_frac(map_2d(1,ainfo%land_index(l))), r_um)
      ! c saturated soil fraction
      c_fsat_soilt(l,1) = real(c_sat_frac(map_2d(1,ainfo%land_index(l))), r_um)
      ! a soil wetness fraction
      a_fwet_soilt(l,1) = real(a_wet_frac(map_2d(1,ainfo%land_index(l))), r_um)
      ! c soil wetness fraction
      c_fwet_soilt(l,1) = real(c_wet_frac(map_2d(1,ainfo%land_index(l))), r_um)
      ! Net primary productivity diagnostic
      npp_gb(l) = real(net_prim_prod(map_2d(1,ainfo%land_index(l))), r_um)
    end do

    ! Prognostics:
    ! Land tile temperatures
    do l = 1, land_pts
      do n = 1, nsurft
        progs%tstar_surft(l, n) = real(tile_temperature(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      end do
    end do

    ! Soil ancillaries and prognostics
    do l = 1, land_pts
      psparms%satcon_soilt(l, 1, 0) = real(soil_cond_sat(map_2d(1,ainfo%land_index(l))), r_um)
      do m = 1, sm_levels
        ! Volumetric soil moisture at wilting point (smvcwt_soilt)
        psparms%smvcwt_soilt(l,1,m) = real(soil_moist_wilt(map_2d(1,ainfo%land_index(l))), r_um)
        ! Volumetric soil moisture at critical point (smvccl_soilt)
        psparms%smvccl_soilt(l,1,m) = real(soil_moist_crit(map_2d(1,ainfo%land_index(l))), r_um)
        ! Volumetric soil moisture at saturation (smvcst_soilt)
        psparms%smvcst_soilt(l,1,m) = real(soil_moist_sat(map_2d(1,ainfo%land_index(l))), r_um)
        ! Saturated soil conductivity (satcon_soilt)
        psparms%satcon_soilt(l,1,m) = real(soil_cond_sat(map_2d(1,ainfo%land_index(l))), r_um)
        ! Soil thermal capacity (hcap_soilt)
        psparms%hcap_soilt(l,1,m) = real(soil_thermal_cap(map_2d(1,ainfo%land_index(l))), r_um)
        ! Saturated soil water suction (sathh_soilt)
        psparms%sathh_soilt(l,1,m) = real(soil_suction_sat(map_2d(1,ainfo%land_index(l))), r_um)
        ! Clapp and Hornberger b coefficient (bexp_soilt)
        psparms%bexp_soilt(l,1,m) = real(clapp_horn_b(map_2d(1,ainfo%land_index(l))), r_um)
        ! Soil temperature (t_soil_soilt)
        progs%t_soil_soilt(l,1,m) = real(soil_temperature(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        ! Soil moisture content (kg m-2, soil_layer_moisture/smcl_soilt)
        progs%smcl_soilt(l,1,m) = real(soil_moisture(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        ! Unfrozen soil moisture proportion (sthu_soilt)
        psparms%sthu_soilt(l,1,m) = real(unfrozen_soil_moisture(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        ! Frozen soil moisture proportion (sthf_soilt)
        psparms%sthf_soilt(l,1,m) = real(frozen_soil_moisture(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        ! Water extraction (ext_soilt [=ext in bl_kernel))
        fluxes%ext_soilt(l,1,m) = real(water_extraction(map_soil(1,ainfo%land_index(l))+m-1), r_um)
      end do
    end do

    ! Soil thermal conductivity (hcon_soilt, hcons_soilt)
    allocate(hcons_soilt(land_pts,nsoilt))
    do l = 1, land_pts
      psparms%hcon_soilt(l,1,:) = real(soil_thermal_cond(map_2d(1,ainfo%land_index(l))), r_um)
      hcons_soilt(l,1) = real(thermal_cond_wet_soil(map_2d(1,ainfo%land_index(l))), r_um)
    end do

    ! Soil and land ice ancils dependant on smvcst_soilt
    ! (soil moisture saturation limit)
    soil_pts = 0
    lice_pts = 0
    do l = 1, land_pts
      if ( psparms%smvcst_soilt(l, 1, 1) > 0.0_r_um ) then
        soil_pts = soil_pts + 1
        ainfo%soil_index(soil_pts) = l
        ainfo%l_soil_point(l) = .true.
      else
        lice_pts = lice_pts + 1
        ainfo%lice_index(lice_pts) = l
        ainfo%l_lice_point(l) = .true.
      end if
    end do

    ! Calculate the infiltration rate
    call infiltration_rate(land_pts, ntiles, ainfo%surft_pts,                  &
                           ainfo%surft_index, psparms%satcon_soilt,            &
                           ainfo%frac_surft, psparms%infil_surft)

    ! Canopy water on each tile (canopy_surft)
    do l = 1, land_pts
      do n = 1, nsurft
        progs%canopy_surft(l, n) = real(canopy_water(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      end do
    end do

    ! Snow prognostics
    do l = 1, land_pts
      do n = 1, nsurft
        ! Lying snow mass on land tiles
        progs%snow_surft(l,n) = real(tile_snow_mass(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        ! Snow grain size on tiles (microns)
        progs%rgrain_surft(l,n) = real(tile_snow_rgrain(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        ! Number of snow layers on tiles (nsnow_surft)
        progs%nsnow_surft(l,n) = n_snow_layers(map_tile(1,ainfo%land_index(l))+n-1)
        ! Snow depth on tiles (snowdepth_surft)
        progs%snowdepth_surft(l,n) = real(snow_depth(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        ! Snow mass under canopy
        progs%snow_grnd_surft(l,n) = real(snow_under_canopy(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        ! Snowpack density (rho_snow_grnd_surft)
        progs%rho_snow_grnd_surft(l,n) = real(snowpack_density(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        do j = 1, nsmax
          i_snow = nsmax*(n-1) + j - 1
          ! Thickness of snow layers
          progs%ds_surft(l,n,j) = real(snow_layer_thickness(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          ! Mass of ice in snow layers
          progs%sice_surft(l,n,j) = real(snow_layer_ice_mass(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          ! Mass of liquid in snow layers
          progs%sliq_surft(l,n,j) = real(snow_layer_liq_mass(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          ! Temperature of snow layers
          progs%tsnow_surft(l,n,j) = real(snow_layer_temp(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          ! Grain size of snow layers
          progs%rgrainl_surft(l,n,j) = real(snow_layer_rgrain(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
        end do
      end do
    end do

    ! Snow melt
    do l = 1, land_pts
      do n = 1, nsurft
        fluxes%melt_surft(l,n) = real(snowice_melt(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      end do
    end do

    allocate(fsat_soilt(land_pts, nsoilt))
    allocate(zw_soilt(land_pts, nsoilt))
    allocate(sthzw_soilt(land_pts, nsoilt))
    do l = 1, land_pts
      ! Soil saturated fraction
      fsat_soilt(l,1) = real(soil_sat_frac(map_2d(1,ainfo%land_index(l))), r_um)
      ! Water table depth
      zw_soilt(l,1) = real(water_table(map_2d(1,ainfo%land_index(l))), r_um)
      ! Soil wetness below soil column
      sthzw_soilt(l,1) = real(wetness_under_soil(map_2d(1,ainfo%land_index(l))), r_um)
    end do

  !----------------------------------------------------------------------------
  ! Call to surf_couple_extra using JULESvn5.4 standalone variable names

    !fqw_surft is not set so will not pick up value calculated via bl_imp
    !However, this is only presently used by river routing so has no effect in
    !LFRic

    !sw_surft is not set in this kernel so will not pick up appropriate values
    !This will affect water resource, irrigation and lakes once coupled to LFRic

    allocate(u_s_std_surft(land_pts, ntiles))
    allocate(gamtot_soilt(land_pts, nsoilt))
    allocate(ti_sig_soilt(land_pts, nsoilt))
    allocate(fwetl_soilt(land_pts, nsoilt))
    allocate(frac_agr_gb(land_pts))
    allocate(dhf_surf_minus_soil(land_pts))
    allocate(tot_surf_runoff(land_pts))
    allocate(tot_sub_runoff(land_pts))
    allocate(inlandout_atm_gb(land_pts))

    call surf_couple_extra(                                                   &
    !Driving data and associated INTENT(IN)
    u_1_ij, v_1_ij,                                                           &

    !Misc INTENT(IN)
    a_step, smlt, ainfo%frac_surft, hcons_soilt, rho_star,                    &

    !IN
    land_pts, seg_len, 1, river_row_length, river_rows,                       &
    ls_graup_ij, cca_2d_ij, nsurft, ainfo%surft_pts,                          &
    lice_pts, soil_pts, stf_sub_surf_roff, fexp_soilt,                        &
    gamtot_soilt, ti_mean_soilt, ti_sig_soilt, flash_rate_ancil,              &
    pop_den_ancil, wealth_index_ancil, a_fsat_soilt, c_fsat_soilt,            &
    a_fwet_soilt, c_fwet_soilt, ntype, delta_lambda, delta_phi,               &
    cos_theta_latitude, aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa,   &
    yua, yva, g_p_field, g_r_field, nproc, global_row_length, global_rows,    &
    global_river_row_length, global_river_rows, flandg, trivdir, trivseq,     &
    r_area, slope, flowobs1, r_inext, r_jnext, r_land, frac_agr_gb,           &
    soil_clay_ij, npp_gb,  u_s_std_surft,                                     &

    !IN OUT
    a_steps_since_riv,                                                        &
    fsat_soilt, fwetl_soilt, zw_soilt, sthzw_soilt,                           &
    ls_rainfrac_gb, substore, surfstore, flowin, bflowin,                     &
    tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,                 &
    asteps_since_triffid,                                                     &
    inlandout_atm_gb,                                                         &

    !OUT- mostly for SCM diagnostics output below
    dhf_surf_minus_soil,                                                      &

    ! IN
    land_sea_mask,                                                            &

    ! JULES TYPES containing field data
    crop_vars, psparms, toppdm, fire_vars, ainfo, trif_vars, soilecosse,      &
    urban_param, progs, trifctltype, coast, jules_vars,                       &
    fluxes,                                                                   &
    lake_vars,                                                                &
    forcing,                                                                  &
    rivers,                                                                   &
    chemvars, water_resources,                                                &
    wtrac_jls,                                                                &
    work_vars_cbl                                                             &
    )

    deallocate(inlandout_atm_gb)
    deallocate(tot_sub_runoff)
    deallocate(tot_surf_runoff)
    deallocate(dhf_surf_minus_soil)
    deallocate(frac_agr_gb)
    deallocate(fwetl_soilt)
    deallocate(ti_sig_soilt)
    deallocate(gamtot_soilt)
    deallocate(u_s_std_surft)

  !---------------------------------------------------------------------------
  ! Return the updated prognostic values to jules_prognostics.

    do n = 1, nsurft
      do l = 1, land_pts
        ! Canopy water on each tile (canopy_surft)
        canopy_water(map_tile(1,ainfo%land_index(l))+n-1) = real(progs%canopy_surft(l,n), r_def)
        ! Lying snow mass on land tiles
        tile_snow_mass(map_tile(1,ainfo%land_index(l))+n-1) = real(progs%snow_surft(l,n), r_def)
        ! Number of snow layers on tiles (nsnow_surft)
        n_snow_layers(map_tile(1,ainfo%land_index(l))+n-1) = progs%nsnow_surft(l,n)
        ! Snow depth on tiles (snowdepth_surft)
        snow_depth(map_tile(1,ainfo%land_index(l))+n-1) = real(progs%snowdepth_surft(l,n), r_def)
        ! Snow grain size on tiles (microns)
        tile_snow_rgrain(map_tile(1,ainfo%land_index(l))+n-1) = real(progs%rgrain_surft(l,n), r_def)
        ! Snow mass under canopy
        snow_under_canopy(map_tile(1,ainfo%land_index(l))+n-1) = real(progs%snow_grnd_surft(l,n), r_def)
        ! Snowpack density (rho_snow_grnd_surft)
        snowpack_density(map_tile(1,ainfo%land_index(l))+n-1) = real(progs%rho_snow_grnd_surft(l,n), r_def)
        ! Total snow and ice melt
        snowice_melt(map_tile(1,ainfo%land_index(l))+n-1) = real(fluxes%melt_surft(l,n), r_def)
      end do
    end do
    do n = 1, nsurft
      do j = 1, nsmax
        do l = 1, land_pts
          i_snow = nsmax*(n-1) + j - 1 
          ! Thickness of snow layers
          snow_layer_thickness(map_snow(1,ainfo%land_index(l))+i_snow) = real(progs%ds_surft(l,n,j), r_def)
          ! Mass of ice in snow layers
          snow_layer_ice_mass(map_snow(1,ainfo%land_index(l))+i_snow) = real(progs%sice_surft(l,n,j), r_def)
          ! Mass of liquid in snow layers
          snow_layer_liq_mass(map_snow(1,ainfo%land_index(l))+i_snow) = real(progs%sliq_surft(l,n,j), r_def)
          ! Temperature of snow layers
          snow_layer_temp(map_snow(1,ainfo%land_index(l))+i_snow) = real(progs%tsnow_surft(l,n,j), r_def)
          ! Grain size of snow layers
          snow_layer_rgrain(map_snow(1,ainfo%land_index(l))+i_snow) = real(progs%rgrainl_surft(l,n,j), r_def)
        end do
      end do
    end do

    do m = 1, sm_levels
      do l = 1, land_pts
        ! Soil temperature (t_soil_soilt)
        soil_temperature(map_soil(1,ainfo%land_index(l))+m-1) = real(progs%t_soil_soilt(l,1,m), r_def)
        ! Soil moisture content (kg m-2, soil_layer_moisture)
        soil_moisture(map_soil(1,ainfo%land_index(l))+m-1) = real(progs%smcl_soilt(l,1,m), r_def)
        ! Unfrozen soil moisture proportion (sthu_soilt)
        unfrozen_soil_moisture(map_soil(1,ainfo%land_index(l))+m-1) = real(psparms%sthu_soilt(l,1,m), r_def)
        ! Frozen soil moisture proportion (sthf_soilt)
        frozen_soil_moisture(map_soil(1,ainfo%land_index(l))+m-1) = real(psparms%sthf_soilt(l,1,m), r_def)
      end do
    end do

    do l = 1, land_pts
      ! Soil saturated fraction
      soil_sat_frac(map_2d(1,ainfo%land_index(l))) = real(fsat_soilt(l,1), r_def)
      ! Water table depth
      water_table(map_2d(1,ainfo%land_index(l))) = real(zw_soilt(l,1), r_def)
      ! Wetness below soil column
      wetness_under_soil(map_2d(1,ainfo%land_index(l))) = real(sthzw_soilt(l,1), r_def)
      ! River runoffs
      surface_runoff(map_2d(1,ainfo%land_index(l))) = real(fluxes%surf_roff_gb(l), r_def)
      sub_surface_runoff(map_2d(1,ainfo%land_index(l))) = real(fluxes%sub_surf_roff_gb(l), r_def)
    end do

    if (.not. associated(soil_moisture_content, empty_real_data) ) then
      do l = 1, land_pts
        soil_moisture_content(map_2d(1,ainfo%land_index(l))) = progs%smc_soilt(l,1)
      end do
    end if

    if (.not. associated(grid_snow_mass, empty_real_data) .and. &
         l_hydrology) then
      do i = 1, seg_len
        grid_snow_mass(map_2d(1,i)) = progs%snow_mass_ij(i,1)
      end do
    end if

    if (.not. associated(throughfall, empty_real_data) ) then
      do n = 1, nsurft
        do l = 1, land_pts
          throughfall(map_tile(1,ainfo%land_index(l))+n-1) = fluxes%tot_tfall_surft(l,n)
        end do
      end do
    end if

    deallocate(ls_rainfrac_gb)
    deallocate(fexp_soilt)
    deallocate(ti_mean_soilt)
    deallocate(a_fsat_soilt)
    deallocate(c_fsat_soilt)
    deallocate(a_fwet_soilt)
    deallocate(c_fwet_soilt)
    deallocate(npp_gb)
    deallocate(hcons_soilt)
    deallocate(fsat_soilt)
    deallocate(zw_soilt)
    deallocate(sthzw_soilt)

    call ancil_info_nullify(ainfo)
    call ancil_info_dealloc(ainfo_data)

    call forcing_nullify(forcing)
    call forcing_dealloc(forcing_data)

    call crop_vars_nullify(crop_vars)
    call crop_vars_dealloc(crop_vars_data)

    call top_pdm_nullify(toppdm)
    call top_pdm_dealloc(top_pdm_data)

    call soil_ecosse_vars_nullify(soilecosse)
    call soil_ecosse_vars_dealloc(soil_ecosse_vars_data)

    call lake_nullify(lake_vars)
    call lake_dealloc(lake_data)

    call fire_vars_nullify(fire_vars)
    call fire_vars_dealloc(fire_vars_data)

    call trifctl_nullify(trifctltype)
    call trifctl_dealloc(trifctl_data)

    call urban_param_nullify(urban_param)
    call urban_param_dealloc(urban_param_data)

    call trif_vars_nullify(trif_vars)
    call trif_vars_dealloc(trif_vars_data)

    call psparms_nullify(psparms)
    call psparms_dealloc(psparms_data)

    call jules_vars_dealloc(jules_vars_data)
    call jules_vars_nullify(jules_vars)

    call prognostics_nullify(progs)
    call prognostics_dealloc(progs_data)

    call fluxes_nullify(fluxes)
    call fluxes_dealloc(fluxes_data)

    call rivers_nullify(rivers)
    call rivers_dealloc(rivers_data)

    call chemvars_nullify(chemvars)
    call chemvars_dealloc(chemvars_data)

    call water_resources_nullify(water_resources)
    call water_resources_dealloc(water_resources_data)

    call wtrac_jls_nullify(wtrac_jls)
    call wtrac_jls_dealloc(wtrac_jls_data)

  end subroutine jules_extra_code

end module jules_extra_kernel_mod
