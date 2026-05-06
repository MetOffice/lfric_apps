!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the explicit JULES surface exchange scheme.
module jules_exp_kernel_mod

  use argument_mod,           only : arg_type,                   &
                                     GH_FIELD, GH_SCALAR,        &
                                     GH_REAL, GH_INTEGER,        &
                                     GH_READ, GH_WRITE,          &
                                     GH_READWRITE, DOMAIN,       &
                                     ANY_DISCONTINUOUS_SPACE_1,  &
                                     ANY_DISCONTINUOUS_SPACE_2,  &
                                     ANY_DISCONTINUOUS_SPACE_3,  &
                                     ANY_DISCONTINUOUS_SPACE_4,  &
                                     ANY_DISCONTINUOUS_SPACE_5,  &
                                     ANY_DISCONTINUOUS_SPACE_6,  &
                                     ANY_DISCONTINUOUS_SPACE_7,  &
                                     ANY_DISCONTINUOUS_SPACE_8,  &
                                     ANY_DISCONTINUOUS_SPACE_9,  &
                                     ANY_DISCONTINUOUS_SPACE_10, &
                                     STENCIL, REGION
  use constants_mod,          only : i_def, i_um, r_def, r_um, rmdi
  use empty_data_mod,         only : empty_real_data
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type
  use radiation_config_mod,   only : topography, topography_horizon
  use jules_radiation_config_mod,  only : l_albedo_obs
  use jules_sea_seaice_config_mod, only : iseasurfalg,                         &
                                     iseasurfalg_specified_roughness,          &
                                     buddy_sea, buddy_sea_on,                  &
                                     z0m_specified_nml => z0m_specified,       &
                                     z0h_specified_nml => z0h_specified
  use water_constants_mod,    only : tfs

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: jules_exp_kernel_type
    private
    type(arg_type) :: meta_args(108) = (/                                      &
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! theta_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3, STENCIL(REGION)),      &! u_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3, STENCIL(REGION)),      &! v_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_v_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_cl_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_cf_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! z0msea_2d
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! z0m_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2, STENCIL(REGION)),&! tile_fraction
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! leaf_area_index
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! canopy_height
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! peak_to_trough_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! silhouette_area_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_albedo
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_roughness
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_wilt
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_crit
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_sat
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_thermal_cond
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_suction_sat
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! clapp_horn_b
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_respiration
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! thermal_cond_wet_soil
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1, STENCIL(REGION)),&! sea_u_current
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1, STENCIL(REGION)),&! sea_v_current
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_conductivity
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_pensolar
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_pensolar_frac_direct
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_pensolar_frac_diffuse
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_snow_mass
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,    ANY_DISCONTINUOUS_SPACE_2),&! n_snow_layers
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! snow_depth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_thickness
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_ice_mass
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_liq_mass
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_temp
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! surface_conductance
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canopy_water
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! soil_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! soil_moisture
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! unfrozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! frozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! tile_heat_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! tile_moisture_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! net_prim_prod
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cos_zen_angle
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! skyview
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! sw_up_tile
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_lw_grey_albedo
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lw_down_surf
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_blue_surf
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_direct_blue_surf
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! dd_mf_cb
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! ozone
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_liquid
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! rhokm_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7),&! surf_interp
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! rhokh_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! moist_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! heat_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! gradrinr
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! alpha1_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! ashtf_prime_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! dtstar_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! fracaero_t_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! fracaero_s_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! z0h_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! z0m_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! rhokh_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! chr1p5m_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! resfs_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! gc_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! canhc_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_8),&! tile_water_extract
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_1),&! blend_height_tq
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! z0m_eff
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! ustar
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_avail
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_3),&! snow_unload_rate
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_9),&! albedo_obs_scaling
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_clay_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_sand_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_10),&! dust_div_mrel
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_10),&! dust_div_flux
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ                             ), &! day_of_year
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ                             ), &! second_of_day
         arg_type(GH_SCALAR, GH_REAL,    GH_READ                             ), &! flux_e
         arg_type(GH_SCALAR, GH_REAL,    GH_READ                             ), &! flux_h
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbwrr
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbhwr
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbhgt
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbztm
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbdisp
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! rhostar
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! recip_l_mo_sea
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! t1_sd
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! q1_sd
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! diag__gross_prim_prod
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! diag__z0h_eff
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,    ANY_DISCONTINUOUS_SPACE_1) &! ocn_cpl_point
         /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: jules_exp_code
  end type

  public :: jules_exp_code

contains

  !> @brief Interface to the Jules explicit scheme
  !> @details The Jules explicit surface exchange does surface forcing of:
  !>             heat, momentum and moisture, and dust emission
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     theta_in_wth           Potential temperature field
  !> @param[in]     exner_in_wth           Exner pressure field in wth space
  !> @param[in]     u_in_w3                'Zonal' wind in density space
  !> @param[in]     v_in_w3                'Meridional' wind in density space
  !> @param[in]     m_v_n                  Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n                 Cloud liquid mixing ratio at time level n
  !> @param[in]     m_cf_n                 Cloud frozen mixing ratio at time level n
  !> @param[in]     height_w3              Height of density space above surface
  !> @param[in]     height_wth             Height of theta space above surface
  !> @param[in]     zh_2d                  Boundary layer depth
  !> @param[in,out] z0msea_2d              Roughness length (sea)
  !> @param[in,out] z0m_2d                 Cell roughness length
  !> @param[in]     tile_fraction          Surface tile fractions
  !> @param[in]     leaf_area_index        Leaf Area Index
  !> @param[in]     canopy_height          Canopy height
  !> @param[in]     peak_to_trough_orog    Half of peak-to-trough height over root(2) of orography
  !> @param[in]     silhouette_area_orog   Silhouette area of orography
  !> @param[in]     soil_albedo            Snow-free soil albedo
  !> @param[in]     soil_roughness         Bare soil surface roughness length
  !> @param[in]     soil_moist_wilt        Volumetric soil moisture at wilting point
  !> @param[in]     soil_moist_crit        Volumetric soil moisture at critical point
  !> @param[in]     soil_moist_sat         Volumetric soil moisture at saturation
  !> @param[in]     soil_thermal_cond      Soil thermal conductivity
  !> @param[in]     soil_suction_sat       Saturated soil water suction
  !> @param[in]     clapp_horn_b           Clapp and Hornberger b coefficient
  !> @param[in,out] soil_respiration       Soil respiration  (kg m-2 s-1)
  !> @param[in,out] thermal_cond_wet_soil  Thermal conductivity of wet soil (W m-1 K-1)
  !> @param[in]     sea_u_current          Sea surface current U component (m/s)
  !> @param[in]     sea_v_current          Sea surface current V component (m/s)
  !> @param[in]     sea_ice_temperature    Bulk temperature of sea-ice (K)
  !> @param[in]     sea_ice_conductivity   Sea ice thermal conductivity (W m-2 K-1)
  !> @param[in,out] sea_ice_pensolar       Sea ice penetrating solar radiation (W m-2)
  !> @param[in] sea_ice_pensolar_frac_direct   Fraction of penetrating solar (direct visible) into sea ice
  !> @param[in] sea_ice_pensolar_frac_diffuse  Fraction of penetrating solar (diffuse visible) into sea ice
  !> @param[in,out] tile_temperature       Surface tile temperatures
  !> @param[in]     tile_snow_mass         Snow mass on tiles (kg/m2)
  !> @param[in]     n_snow_layers          Number of snow layers on tiles
  !> @param[in]     snow_depth             Snow depth on tiles
  !> @param[in]     snow_layer_thickness   Thickness of snow layers (m)
  !> @param[in]     snow_layer_ice_mass    Mass of ice in snow layers (kg m-2)
  !> @param[in]     snow_layer_liq_mass    Mass of liquid in snow layers (kg m-2)
  !> @param[in]     snow_layer_temp        Temperature of snow layer (K)
  !> @param[in,out] surface_conductance    Surface conductance
  !> @param[in]     canopy_water           Canopy water on each tile
  !> @param[in]     soil_temperature       Soil temperature
  !> @param[in]     soil_moisture          Soil moisture content (kg m-2)
  !> @param[in]     unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in]     frozen_soil_moisture   Frozen soil moisture proportion
  !> @param[in,out] tile_heat_flux         Surface heat flux
  !> @param[in,out] tile_moisture_flux     Surface moisture flux
  !> @param[in,out] net_prim_prod          Net Primary Productivity
  !> @param[in]     cos_zen_angle          Cosine of solar zenith angle
  !> @param[in]     skyview                Skyview / area enhancement factor
  !> @param[in]     sw_up_tile             Upwelling SW radiation on surface tiles
  !> @param[in]     tile_lw_grey_albedo    Surface tile longwave grey albedo
  !> @param[in]     sw_down_surf           Downwelling SW radiation at surface
  !> @param[in]     lw_down_surf           Downwelling LW radiation at surface
  !> @param[in]     sw_down_blue_surf      Photosynthetically active SW down
  !> @param[in]     sw_direct_blue_surf    Downwelling direct visible light
  !> @param[in]     dd_mf_cb               Downdraft massflux at cloud base (Pa/s)
  !> @param[in]     ozone                  Ozone field
  !> @param[in]     cf_bulk                Bulk cloud fraction
  !> @param[in]     cf_liquid              Liquid cloud fraction
  !> @param[in,out] rhokm_bl               Momentum eddy diffusivity on BL levels
  !> @param[in,out] surf_interp            Surface variables for regridding
  !> @param[in,out] rhokh_bl               Heat eddy diffusivity on BL levels
  !> @param[in,out] moist_flux_bl          Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl           Vertical heat flux on BL levels
  !> @param[in,out] gradrinr               Gradient Richardson number in wth
  !> @param[in,out] alpha1_tile            dqsat/dT in surface layer on tiles
  !> @param[in,out] ashtf_prime_tile       Heat flux coefficient on tiles
  !> @param[in,out] dtstar_tile            Change in surface temperature on tiles
  !> @param[in,out] fracaero_t_tile        Total fraction of moisture flux with only aerodynamic resistance
  !> @param[in,out] fracaero_s_tile        Fraction of moisture flux with only aerodynamic resistance from frozen surface
  !> @param[in,out] z0h_tile               Heat roughness length on tiles
  !> @param[in,out] z0m_tile               Momentum roughness length on tiles
  !> @param[in,out] rhokh_tile             Surface heat diffusivity on tiles
  !> @param[in,out] chr1p5m_tile           1.5m transfer coefficients on tiles
  !> @param[in,out] resfs_tile             Combined aerodynamic resistance
  !> @param[in,out] gc_tile                Stomatal conductance on tiles (m s-1)
  !> @param[in,out] canhc_tile             Canopy heat capacity on tiles
  !> @param[in,out] tile_water_extract     Extraction of water from each tile
  !> @param[in,out] blend_height_tq        Blending height for wth levels
  !> @param[in,out] z0m_eff                Grid mean effective roughness length
  !> @param[in,out] ustar                  Friction velocity
  !> @param[in,out] soil_moist_avail       Available soil moisture for evaporation
  !> @param[in,out] snow_unload_rate       Unloading of snow from PFTs by wind
  !> @param[in]     albedo_obs_scaling     Scaling factor to adjust albedos by
  !> @param[in]     soil_clay_2d           Soil clay fraction
  !> @param[in]     soil_sand_2d           Soil sand fraction
  !> @param[in]     dust_div_mrel          Relative soil mass in CLASSIC size divisions
  !> @param[in,out] dust_div_flux          Dust emission fluxes in CLASSIC size divisions (kg m-2 s-1)
  !> @param[in]     day_of_year            The day of the year
  !> @param[in]     second_of_day          The second of the day
  !> @param[in]     flux_e                 Latent heat flux
  !> @param[in]     flux_h                 Sensible heat flux
  !> @param[in]     urbwrr                 Urban repeating width ratio
  !> @param[in]     urbhwr                 Urban height-to-width ratio
  !> @param[in]     urbhgt                 Urban building height
  !> @param[in]     urbztm                 Urban effective roughness length
  !> @param[in]     urbdisp                Urban displacement height
  !> @param[in,out] rhostar_2d             Surface density
  !> @param[in,out] recip_l_mo_sea_2d      Inverse Obukhov length over sea only
  !> @param[in,out] t1_sd_2d               StDev of level 1 temperature
  !> @param[in,out] q1_sd_2d               StDev of level 1 humidity
  !> @param[in,out] gross_prim_prod        Diagnostic: Gross Primary Productivity
  !> @param[in,out] z0h_eff                Diagnostic: Gridbox mean effective roughness length for scalars
  !> @param[in,out] ocn_cpl_point          Diagnostic: Coupling point mask
  !> @param[in]     ndf_wth                Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth               Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth                Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3                 Number of DOFs per cell for density space
  !> @param[in]     undf_w3                Number of unique DOFs for density space
  !> @param[in]     map_w3                 Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d                 Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d                Number of unique DOFs for 2D fields
  !> @param[in]     map_2d                 Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile               Number of DOFs per cell for tiles
  !> @param[in]     undf_tile              Number of total DOFs for tiles
  !> @param[in]     map_tile               Dofmap for cell for surface tiles
  !> @param[in]     ndf_pft                Number of DOFs per cell for PFTs
  !> @param[in]     undf_pft               Number of total DOFs for PFTs
  !> @param[in]     map_pft                Dofmap for cell for PFTs
  !> @param[in]     ndf_sice               Number of DOFs per cell for sice
  !> @param[in]     undf_sice              Number of total DOFs for sice
  !> @param[in]     map_sice               Dofmap for cell for sice
  !> @param[in]     ndf_snow               Number of DOFs per cell for snow
  !> @param[in]     undf_snow              Number of total DOFs for snow
  !> @param[in]     map_snow               Dofmap for cell for snow
  !> @param[in]     ndf_soil               Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil              Number of total DOFs for soil levels
  !> @param[in]     map_soil               Dofmap for cell for soil levels
  !> @param[in]     ndf_surf               Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf              Number of unique DOFs for surface variables
  !> @param[in]     map_surf               Dofmap for the cell at the base of the column for surface variables
  !> @param[in]     ndf_smtile             Number of DOFs per cell for soil levels and tiles
  !> @param[in]     undf_smtile            Number of total DOFs for soil levels and tiles
  !> @param[in]     map_smtile             Dofmap for cell for soil levels and tiles
  !> @param[in]     ndf_scal               Number of DOFs per cell for albedo scaling
  !> @param[in]     undf_scal              Number of total DOFs for albedo scaling
  !> @param[in]     map_scal               Dofmap for cell at the base of the column
  !> @param[in]     ndf_dust               Number of DOFs per cell for dust divisions
  !> @param[in]     undf_dust              Number of total DOFs for dust divisions
  !> @param[in]     map_dust               Dofmap for cell for dust divisions
  subroutine jules_exp_code(nlayers, seg_len, seg_len_halo,       &
                           theta_in_wth,                          &
                           exner_in_wth,                          &
                           u_in_w3,                               &
                           u_w3_stencil_size, u_w3_stencil,       &
                           v_in_w3,                               &
                           v_w3_stencil_size, v_w3_stencil,       &
                           m_v_n,                                 &
                           m_cl_n,                                &
                           m_cf_n,                                &
                           height_w3,                             &
                           height_wth,                            &
                           zh_2d,                                 &
                           z0msea_2d,                             &
                           z0m_2d,                                &
                           tile_fraction,                         &
                           tile_stencil_size, tile_stencil,       &
                           leaf_area_index,                       &
                           canopy_height,                         &
                           peak_to_trough_orog,                   &
                           silhouette_area_orog,                  &
                           soil_albedo,                           &
                           soil_roughness,                        &
                           soil_moist_wilt,                       &
                           soil_moist_crit,                       &
                           soil_moist_sat,                        &
                           soil_thermal_cond,                     &
                           soil_suction_sat,                      &
                           clapp_horn_b,                          &
                           soil_respiration,                      &
                           thermal_cond_wet_soil,                 &
                           sea_u_current, sea_u_w3_stencil_size,  &
                           sea_u_w3_stencil,                      &
                           sea_v_current, sea_v_w3_stencil_size,  &
                           sea_v_w3_stencil,                      &
                           sea_ice_temperature,                   &
                           sea_ice_conductivity,                  &
                           sea_ice_pensolar,                      &
                           sea_ice_pensolar_frac_direct,          &
                           sea_ice_pensolar_frac_diffuse,         &
                           tile_temperature,                      &
                           tile_snow_mass,                        &
                           n_snow_layers,                         &
                           snow_depth,                            &
                           snow_layer_thickness,                  &
                           snow_layer_ice_mass,                   &
                           snow_layer_liq_mass,                   &
                           snow_layer_temp,                       &
                           surface_conductance,                   &
                           canopy_water,                          &
                           soil_temperature,                      &
                           soil_moisture,                         &
                           unfrozen_soil_moisture,                &
                           frozen_soil_moisture,                  &
                           tile_heat_flux,                        &
                           tile_moisture_flux,                    &
                           net_prim_prod,                         &
                           cos_zen_angle,                         &
                           skyview,                               &
                           sw_up_tile,                            &
                           tile_lw_grey_albedo,                   &
                           sw_down_surf,                          &
                           lw_down_surf,                          &
                           sw_down_blue_surf,                     &
                           sw_direct_blue_surf,                   &
                           dd_mf_cb,                              &
                           ozone,                                 &
                           cf_bulk,                               &
                           cf_liquid,                             &
                           rhokm_bl,                              &
                           surf_interp,                           &
                           rhokh_bl,                              &
                           moist_flux_bl,                         &
                           heat_flux_bl,                          &
                           gradrinr,                              &
                           alpha1_tile,                           &
                           ashtf_prime_tile,                      &
                           dtstar_tile,                           &
                           fracaero_t_tile,                       &
                           fracaero_s_tile,                       &
                           z0h_tile,                              &
                           z0m_tile,                              &
                           rhokh_tile,                            &
                           chr1p5m_tile,                          &
                           resfs_tile,                            &
                           gc_tile,                               &
                           canhc_tile,                            &
                           tile_water_extract,                    &
                           blend_height_tq,                       &
                           z0m_eff,                               &
                           ustar,                                 &
                           soil_moist_avail,                      &
                           snow_unload_rate,                      &
                           albedo_obs_scaling,                    &
                           soil_clay_2d,                          &
                           soil_sand_2d,                          &
                           dust_div_mrel,                         &
                           dust_div_flux,                         &
                           day_of_year,                           &
                           second_of_day,                         &
                           flux_e,                                &
                           flux_h,                                &
                           urbwrr,                                &
                           urbhwr,                                &
                           urbhgt,                                &
                           urbztm,                                &
                           urbdisp,                               &
                           rhostar_2d,                            &
                           recip_l_mo_sea_2d,                     &
                           t1_sd_2d,                              &
                           q1_sd_2d,                              &
                           gross_prim_prod,                       &
                           z0h_eff,                               &
                           ocn_cpl_point,                         &
                           ndf_wth, undf_wth, map_wth,            &
                           ndf_w3, undf_w3, map_w3,               &
                           ndf_2d, undf_2d, map_2d,               &
                           ndf_tile, undf_tile, map_tile,         &
                           ndf_pft, undf_pft, map_pft,            &
                           ndf_sice, undf_sice, map_sice,         &
                           ndf_snow, undf_snow, map_snow,         &
                           ndf_soil, undf_soil, map_soil,         &
                           ndf_surf, undf_surf, map_surf,         &
                           ndf_smtile, undf_smtile, map_smtile,   &
                           ndf_scal, undf_scal, map_scal,         &
                           ndf_dust, undf_dust, map_dust )

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use gas_calc_all_mod, only: co2_mix_ratio_now
    use jules_control_init_mod, only: n_land_tile, n_sea_ice_tile,        &
                                      first_sea_tile, first_sea_ice_tile, &
                                      n_surf_tile

    !---------------------------------------
    ! UM/Jules modules containing switches or global constants
    !---------------------------------------
    use ancil_info, only: rad_nband, nsoilt, dim_cslayer, nmasst
    use atm_fields_bounds_mod, only: pdims_s, pdims
    use atm_step_local, only: dim_cs1, co2_dim_len, co2_dim_row
    use dust_parameters_mod, only: ndiv, ndivh, ndivl, l_dust_flux_only
    use jules_deposition_mod, only: l_deposition
    use jules_irrig_mod, only: irr_crop, irr_crop_doell
    use jules_sea_seaice_mod, only: nice_use, l_ctile, l_sice_swpen
    use jules_snow_mod, only: nsmax, cansnowtile
    use jules_soil_mod, only: ns_deep, l_bedrock
    use jules_soil_biogeochem_mod, only: dim_ch4layer, soil_bgc_model,         &
                                         soil_model_ecosse, l_layeredc
    use jules_surface_mod, only: l_urban2t, l_flake_model
    use jules_surface_types_mod, only: npft, ntype, ncpft, nnpft, soil
    use jules_urban_mod, only: l_moruses
    use jules_vegetation_mod, only: l_crop, l_triffid, l_phenol, l_use_pft_psi,&
                                    can_rad_mod, l_acclim, l_sugar, l_red
    use jules_water_tracers_mod, only: l_wtrac_jls, n_wtrac_jls, n_evap_srce
    use nlsizes_namelist_mod, only: sm_levels, ntiles, bl_levels
    use planet_constants_mod, only: p_zero, kappa, planet_radius, cp, g, grcp, &
                                    c_virtual, repsilon, r, lcrcp, lsrcp, vkman
    use rad_input_mod, only: co2_mmr
    use bl_option_mod, only: one_third, flux_bc_opt,interactive_fluxes,  &
                             specified_fluxes_only, specified_fluxes_cd, &
                             l_noice_in_turb
    use water_constants_mod, only: lc
    use c_elevate, only: l_elev_absolute_height

    ! subroutines used
    use buoy_tq_mod, only: buoy_tq
    use dust_calc_emiss_frac_mod, only: dust_calc_emiss_frac
    use dust_srce_mod, only: dust_srce
    use qsat_mod, only: qsat_mix
    use surf_couple_explicit_mod, only: surf_couple_explicit
    use sf_diags_mod, only: sf_diag, dealloc_sf_expl, alloc_sf_expl
    use sparm_mod, only: sparm
    use tilepts_mod, only: tilepts

    !---------------------------------------
    ! JULES modules
    !---------------------------------------
    use crop_vars_mod,            only: crop_vars_type, crop_vars_data_type,   &
                                        crop_vars_alloc, crop_vars_assoc,      &
                                        crop_vars_dealloc, crop_vars_nullify
    use prognostics,              only: progs_data_type, progs_type,           &
                                        prognostics_alloc, prognostics_assoc,  &
                                        prognostics_dealloc, prognostics_nullify
    use jules_vars_mod,           only: jules_vars_type, jules_vars_data_type, &
                                        jules_vars_alloc, jules_vars_assoc,    &
                                        jules_vars_dealloc, jules_vars_nullify
    use p_s_parms,                only: psparms_type, psparms_data_type,       &
                                        psparms_alloc, psparms_assoc,          &
                                        psparms_dealloc, psparms_nullify
    use trif_vars_mod,            only: trif_vars_type, trif_vars_data_type,   &
                                        trif_vars_assoc, trif_vars_alloc,      &
                                        trif_vars_dealloc, trif_vars_nullify
    use aero,                     only: aero_type, aero_data_type,             &
                                        aero_assoc, aero_alloc,                &
                                        aero_dealloc, aero_nullify
    use urban_param_mod,          only: urban_param_type,                      &
                                        urban_param_data_type,                 &
                                        urban_param_assoc, urban_param_alloc,  &
                                        urban_param_dealloc, urban_param_nullify
    use trifctl,                  only: trifctl_type, trifctl_data_type,       &
                                        trifctl_assoc, trifctl_alloc,          &
                                        trifctl_dealloc, trifctl_nullify
    use coastal,                  only: coastal_type, coastal_data_type,       &
                                        coastal_assoc, coastal_alloc,          &
                                        coastal_dealloc, coastal_nullify
    use lake_mod,                 only: lake_type, lake_data_type,             &
                                        lake_assoc, lake_alloc,                &
                                        lake_dealloc, lake_nullify
    use ancil_info,               only: ainfo_type, ainfo_data_type,           &
                                        ancil_info_assoc, ancil_info_alloc,    &
                                        ancil_info_dealloc, ancil_info_nullify
    use jules_forcing_mod,        only: forcing_type, forcing_data_type,       &
                                        forcing_assoc, forcing_alloc,          &
                                        forcing_dealloc, forcing_nullify
    use fluxes_mod,               only: fluxes_type, fluxes_data_type,         &
                                        fluxes_alloc, fluxes_assoc,            &
                                        fluxes_nullify, fluxes_dealloc
    use jules_chemvars_mod,       only: chemvars_type, chemvars_data_type,     &
                                        chemvars_alloc, chemvars_assoc,        &
                                        chemvars_dealloc, chemvars_nullify
    use jules_wtrac_type_mod,     only: jls_wtrac_type, jls_wtrac_data_type,   &
                                        wtrac_jls_alloc, wtrac_jls_assoc,      &
                                        wtrac_jls_dealloc, wtrac_jls_nullify
    use progs_cbl_vars_mod, only: progs_cbl_vars_type
    use work_vars_mod_cbl, only: work_vars_type

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len, seg_len_halo
    integer(kind=i_def), intent(in) :: day_of_year, second_of_day
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth,seg_len)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3,seg_len)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d,seg_len)

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile,seg_len)
    integer(kind=i_def), intent(in) :: ndf_pft, undf_pft
    integer(kind=i_def), intent(in) :: map_pft(ndf_pft,seg_len)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil,seg_len)
    integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
    integer(kind=i_def), intent(in) :: map_sice(ndf_sice,seg_len)
    integer(kind=i_def), intent(in) :: ndf_snow, undf_snow
    integer(kind=i_def), intent(in) :: map_snow(ndf_snow,seg_len)
    integer(kind=i_def), intent(in) :: ndf_smtile, undf_smtile
    integer(kind=i_def), intent(in) :: map_smtile(ndf_smtile,seg_len)

    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf
    integer(kind=i_def), intent(in) :: map_surf(ndf_surf,seg_len)
    integer(kind=i_def), intent(in) :: ndf_scal, undf_scal
    integer(kind=i_def), intent(in) :: map_scal(ndf_scal,seg_len)
    integer(kind=i_def), intent(in) :: ndf_dust, undf_dust
    integer(kind=i_def), intent(in) :: map_dust(ndf_dust,seg_len)

    integer(kind=i_def), intent(in) :: u_w3_stencil_size(seg_len), v_w3_stencil_size(seg_len)
    integer(kind=i_def), dimension(ndf_w3,maxval(u_w3_stencil_size),seg_len_halo), intent(in) :: u_w3_stencil
    integer(kind=i_def), dimension(ndf_w3,maxval(v_w3_stencil_size),seg_len_halo), intent(in) :: v_w3_stencil

    integer(kind=i_def), intent(in) :: tile_stencil_size(seg_len)
    integer(kind=i_def), dimension(ndf_tile,maxval(tile_stencil_size),seg_len_halo), intent(in) :: tile_stencil

    real(kind=r_def), dimension(undf_wth), intent(inout):: rhokm_bl,           &
                                                           gradrinr
    real(kind=r_def), dimension(undf_w3),  intent(inout):: rhokh_bl,           &
                                                           moist_flux_bl,      &
                                                           heat_flux_bl
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: u_in_w3, v_in_w3,   &
                                                           height_w3
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           exner_in_wth,       &
                                                           m_v_n, m_cl_n,      &
                                                           m_cf_n,             &
                                                           height_wth,         &
                                                           cf_bulk, cf_liquid, &
                                                           ozone

    real(kind=r_def), dimension(undf_2d), intent(in)    :: zh_2d
    real(kind=r_def), dimension(undf_2d), intent(inout) :: z0msea_2d,          &
                                                           z0m_2d,             &
                                                           z0m_eff,            &
                                                           ustar,              &
                                                           soil_moist_avail,   &
                                                           recip_l_mo_sea_2d,  &
                                                           rhostar_2d,         &
                                                           t1_sd_2d, q1_sd_2d
    integer(kind=i_def), dimension(undf_2d), intent(inout) :: blend_height_tq

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(in) :: tile_snow_mass(undf_tile)
    integer(kind=i_def), intent(in) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(in) :: snow_depth(undf_tile)
    real(kind=r_def), intent(in) :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_heat_flux(undf_tile)
    real(kind=r_def), intent(inout) :: tile_moisture_flux(undf_tile)
    real(kind=r_def), intent(in) :: sw_up_tile(undf_tile)
    real(kind=r_def), intent(in) :: tile_lw_grey_albedo(undf_tile)
    real(kind=r_def), intent(in) :: albedo_obs_scaling(undf_scal)

    real(kind=r_def), intent(in) :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in) :: canopy_height(undf_pft)

    real(kind=r_def), intent(in) :: sea_u_current(undf_2d)
    real(kind=r_def), intent(in) :: sea_v_current(undf_2d)
    integer(kind=i_def), intent(in) :: sea_u_w3_stencil_size(seg_len), sea_v_w3_stencil_size(seg_len)
    integer(kind=i_def), dimension(ndf_w3,maxval(sea_u_w3_stencil_size),seg_len_halo), intent(in) :: sea_u_w3_stencil
    integer(kind=i_def), dimension(ndf_w3,maxval(sea_v_w3_stencil_size),seg_len_halo), intent(in) :: sea_v_w3_stencil

    real(kind=r_def), intent(in) :: sea_ice_temperature(undf_sice)

    real(kind=r_def), intent(in) :: sea_ice_conductivity(undf_sice)
    real(kind=r_def), intent(inout) :: sea_ice_pensolar(undf_sice)
    real(kind=r_def), intent(in) :: sea_ice_pensolar_frac_direct(undf_sice)
    real(kind=r_def), intent(in) :: sea_ice_pensolar_frac_diffuse(undf_sice)

    real(kind=r_def), intent(in) :: peak_to_trough_orog(undf_2d)
    real(kind=r_def), intent(in) :: silhouette_area_orog(undf_2d)
    real(kind=r_def), intent(in) :: soil_albedo(undf_2d)
    real(kind=r_def), intent(in) :: soil_roughness(undf_2d)
    real(kind=r_def), intent(in) :: soil_thermal_cond(undf_2d)
    real(kind=r_def), intent(inout) :: surface_conductance(undf_2d)
    real(kind=r_def), intent(in) :: cos_zen_angle(undf_2d)
    real(kind=r_def), intent(in) :: skyview(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: lw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_blue_surf(undf_2d)
    real(kind=r_def), intent(in) :: sw_direct_blue_surf(undf_2d)
    real(kind=r_def), intent(in) :: dd_mf_cb(undf_2d)
    real(kind=r_def), intent(inout) :: net_prim_prod(undf_2d)
    real(kind=r_def), pointer, intent(inout) :: soil_respiration(:)
    real(kind=r_def), intent(inout) :: thermal_cond_wet_soil(undf_2d)
    integer(kind=i_def), intent(in) :: ocn_cpl_point(undf_2d)

    ! Urban morphology fields (surface_fields)
    real(kind=r_def), intent(in) :: urbwrr(undf_2d)
    real(kind=r_def), intent(in) :: urbhwr(undf_2d)
    real(kind=r_def), intent(in) :: urbhgt(undf_2d)
    real(kind=r_def), intent(in) :: urbztm(undf_2d)
    real(kind=r_def), intent(in) :: urbdisp(undf_2d)

    real(kind=r_def), intent(in) :: soil_moist_wilt(undf_2d)
    real(kind=r_def), intent(in) :: soil_moist_crit(undf_2d)
    real(kind=r_def), intent(in) :: soil_moist_sat(undf_2d)
    real(kind=r_def), intent(in) :: soil_suction_sat(undf_2d)
    real(kind=r_def), intent(in) :: clapp_horn_b(undf_2d)
    real(kind=r_def), intent(in) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(in) :: soil_moisture(undf_soil)
    real(kind=r_def), intent(in) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(in) :: frozen_soil_moisture(undf_soil)

    real(kind=r_def), intent(in) :: snow_layer_thickness(undf_snow)
    real(kind=r_def), intent(in) :: snow_layer_ice_mass(undf_snow)
    real(kind=r_def), intent(in) :: snow_layer_liq_mass(undf_snow)
    real(kind=r_def), intent(in) :: snow_layer_temp(undf_snow)

    real(kind=r_def), intent(inout) :: tile_water_extract(undf_smtile)

    real(kind=r_def), dimension(undf_pft),  intent(inout)  :: snow_unload_rate
    real(kind=r_def), dimension(undf_surf), intent(inout)  :: surf_interp
    real(kind=r_def), dimension(undf_tile), intent(inout):: alpha1_tile,      &
                                                            ashtf_prime_tile, &
                                                            dtstar_tile,      &
                                                            fracaero_t_tile,  &
                                                            fracaero_s_tile,  &
                                                            z0h_tile,         &
                                                            z0m_tile,         &
                                                            rhokh_tile,       &
                                                            chr1p5m_tile,     &
                                                            resfs_tile,       &
                                                            gc_tile,          &
                                                            canhc_tile

    real(kind=r_def), dimension(undf_2d),   intent(in)     :: soil_clay_2d
    real(kind=r_def), dimension(undf_2d),   intent(in)     :: soil_sand_2d
    real(kind=r_def), dimension(undf_dust), intent(in)     :: dust_div_mrel
    real(kind=r_def), dimension(undf_dust), intent(inout)  :: dust_div_flux

    real(kind=r_def), pointer, intent(inout) :: z0h_eff(:), gross_prim_prod(:)
    real(kind=r_def), intent(in) :: flux_h
    real(kind=r_def), intent(in) :: flux_e

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    real(kind=r_um), allocatable :: qs_star(:,:)
    real(r_um) :: rholem, tv1_sd, w_m, dqsdt_star, wthvbar, ch, theta1
    integer(i_def) :: k, i, i_tile, i_sice, n, i_snow, j, l, idiv, m, &
         land_field, ssi_pts, sea_pts

    ! local switches and scalars
    logical, parameter :: l_aero_classic=.false.
    logical :: l_spec_z0
    real(r_def) :: sw_diffuse_blue_surf

    ! profile fields from level 0 upwards
    real(r_um), dimension(seg_len,1,0:1) :: p_theta_levels, q, qcl, qcf

    real(r_um), dimension(co2_dim_len,co2_dim_row) :: co2

    ! single level real fields
    real(r_um), dimension(seg_len,1) ::                                      &
         fqw, ftl, rhokh, ddmfx,                                             &
         z0h_specified, z0m_specified, soil_clay, t1_sd, q1_sd, fb_surf,     &
         rib_gb, vshr, ustargbm, photosynth_act_rad, tstar_land, dtstar_sea, &
         tstar_sice, alpha1_sea, ashtf_prime_sea, chr1p5m_sice, rhokh_sea,   &
         z0hssi, z0mssi, z1_uv_top, z1_tq_top, rhostar, recip_l_mo_sea, sky
    real(r_um), dimension(seg_len,1,1) :: temperature, bt_blend, bq_blend,   &
         bulk_cloud_fraction

    real(r_um), dimension(seg_len_halo,1) ::  flandg,               &
         flandfac, fseafac, cdr10m, rhokm_land, rhokm_ssi, rhokm

    ! single level integer fields
    integer(i_um), dimension(seg_len,1) ::  k_blend_tq
    integer(i_um), dimension(seg_len_halo,1) :: k_blend_uv

    ! single level logical fields
    logical, dimension(seg_len,1) :: land_sea_mask

    ! fields on sea-ice categories
    real(r_um), dimension(seg_len,1,nice_use) ::                             &
         alpha1_sice, ashtf_prime, rhokh_sice, dtstar_sice

    ! field on land points and soil levels
    real(r_um), dimension(:,:), allocatable :: soil_layer_moisture

    ! fields on land points
    real(r_um), dimension(:), allocatable ::                                 &
         hcons_soilt, rhostar_land, sand_land, clay_land, emis_soil,         &
         dust_emsc_land

    ! fields on land points and tiles
    real(r_um), dimension(:,:), allocatable ::                               &
         dtstar_surft, alpha1, ashtf_prime_surft, chr1p5m, resfs,            &
         rhokh_surft, canhc_surft, fracaero_t, fracaero_s,                   &
         epot_surft, flake, resft

    ! field on surface tiles and soil levels
    real(r_um), dimension(:,:,:), allocatable :: wt_ext_surft

    ! Dust emission fields
    real(r_um), dimension(seg_len,1,ndiv) :: dust_flux
    real(r_um), dimension(:,:,:), allocatable :: u_s_t_dry_tile,        &
         u_s_t_tile
    real(r_um), dimension(:,:,:), allocatable :: dust_flux_surft
    real(r_um), dimension(:,:), allocatable :: mrel_land, dust_emiss_frac

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    integer(i_um) :: asteps_since_triffid, ndry_dep_species,                   &
                     river_row_length_dum, river_rows_dum

    real(r_um), dimension(seg_len,1,1) ::                                      &
         bt, bq, bt_cld, bq_cld, a_qs, a_dqsdt, dqsdt
    real(r_um), dimension(seg_len,1) :: charnock_w

    ! This is an idealised fixed value for ustar.
    real(r_um), parameter :: ustar_fixed_value = 0.1_r_um

    real(r_um), dimension(seg_len,1,nice_use) :: radnet_sice

    real(r_um), dimension(:), allocatable :: resp_s_tot_soilt

    !-----------------------------------------------------------------------
    ! JULES Types
    !-----------------------------------------------------------------------
    type(crop_vars_type) :: crop_vars
    type(crop_vars_data_type) :: crop_vars_data
    type(progs_type) :: progs
    type(progs_data_type) :: progs_data
    type(psparms_type) :: psparms
    type(psparms_data_type) :: psparms_data
    type(trif_vars_type) :: trif_vars
    type(trif_vars_data_type) :: trif_vars_data
    type(aero_type) :: aerotype
    type(aero_data_type) :: aero_data
    type(urban_param_type) :: urban_param
    type(urban_param_data_type) :: urban_param_data
    type(trifctl_type) :: trifctltype
    type(trifctl_data_type) :: trifctl_data
    type(coastal_type) :: coast
    type(coastal_data_type) :: coastal_data
    type(lake_type) :: lake_vars
    type(lake_data_type) :: lake_data
    type(ainfo_type) :: ainfo
    type(ainfo_data_type) :: ainfo_data
    type(forcing_type) :: forcing
    type(forcing_data_type) :: forcing_data
    type(fluxes_type) :: fluxes
    type(fluxes_data_type) :: fluxes_data
    type(chemvars_type) :: chemvars
    type(chemvars_data_type) :: chemvars_data
    type(jls_wtrac_type) :: wtrac_jls
    type(jls_wtrac_data_type) :: wtrac_jls_data
    type(jules_vars_type) :: jules_vars
    type(jules_vars_data_type), TARGET :: jules_vars_data
    ! Unused types needed for argument list
    type(progs_cbl_vars_type) :: progs_cbl_vars
    type(work_vars_type) :: work_cbl

    !-----------------------------------------------------------------------
    ! Initialisation of JULES data and pointer types
    !-----------------------------------------------------------------------
    ! Land tile fractions
    land_field = 0
    do i = 1, seg_len_halo
      flandg(i,1) = 0.0_r_um
      do n = 1, n_land_tile
        flandg(i,1) = flandg(i,1) + real(tile_fraction(tile_stencil(1,1,i)+n-1), r_um)
      end do
    end do
    do i = 1, seg_len
      if (flandg(i,1) > 0.0_r_um) then
        land_field = land_field + 1
        land_sea_mask(i,1) = .true.
      else
        land_sea_mask(i,1) = .false.
      end if
    end do

    call crop_vars_alloc(land_field, seg_len, 1, n_land_tile, ncpft,nsoilt,    &
                         sm_levels, l_crop, irr_crop, irr_crop_doell,          &
                         crop_vars_data)
    call crop_vars_assoc(crop_vars, crop_vars_data)

    call prognostics_alloc(land_field, seg_len, 1, n_land_tile, npft, nsoilt,  &
                           sm_levels, ns_deep, nsmax, dim_cslayer, dim_cs1,    &
                           dim_ch4layer, nice_use, nice_use, soil_bgc_model,   &
                           soil_model_ecosse, l_layeredc, l_triffid, l_phenol, &
                           l_bedrock, l_red, nmasst, nnpft, l_acclim,          &
                           l_sugar, progs_data)
    call prognostics_assoc(progs,progs_data)

    call psparms_alloc(land_field, seg_len, 1, nsoilt, sm_levels, dim_cslayer, &
                       n_land_tile, npft, soil_bgc_model, soil_model_ecosse,   &
                       l_use_pft_psi, psparms_data)
    call psparms_assoc(psparms, psparms_data)

    call trif_vars_alloc(land_field, npft, dim_cslayer, nsoilt, dim_cs1,       &
                         l_triffid, l_phenol, trif_vars_data)
    call trif_vars_assoc(trif_vars, trif_vars_data)

    call aero_alloc(land_field, seg_len, 1, n_land_tile,ndiv, aero_data)
    call aero_assoc(aerotype, aero_data)

    call urban_param_alloc(land_field, l_urban2t, l_moruses, urban_param_data)
    call urban_param_assoc(urban_param, urban_param_data)

    call trifctl_alloc(land_field,npft,dim_cslayer,dim_cs1,nsoilt,trifctl_data)
    call trifctl_assoc(trifctltype, trifctl_data)

    call coastal_alloc(land_field, seg_len, 1, seg_len, 1, seg_len, 1,         &
                       nice_use,nice_use,coastal_data,stencil_i_in=maxval(tile_stencil_size),stencil_j_in=1)
    call coastal_assoc(coast, coastal_data)

    call lake_alloc(land_field, l_flake_model, lake_data)
    call lake_assoc(lake_vars, lake_data)

    call ancil_info_alloc(land_field, seg_len, 1, nice_use, nsoilt, ntype,     &
                          ainfo_data)
    call ancil_info_assoc(ainfo, ainfo_data)

    call forcing_alloc(seg_len, 1, seg_len, 1, seg_len, 1, forcing_data)
    call forcing_assoc(forcing, forcing_data)

    call fluxes_alloc(land_field, seg_len, 1, n_land_tile, npft, nsoilt,       &
                      sm_levels, nice_use, nice_use, fluxes_data)
    call fluxes_assoc(fluxes, fluxes_data)

    ! Set num dry dep species to fixed = 1 for now
    ndry_dep_species = 1
    call chemvars_alloc(land_field, seg_len, 1, npft, ntype,                   &
                        l_deposition, ndry_dep_species, chemvars_data)
    call chemvars_assoc(chemvars, chemvars_data)

    ! Set river size to 1 here as the fields with these dimensions are not
    ! used here
    river_row_length_dum = 1
    river_rows_dum = 1
    call wtrac_jls_alloc(land_field, seg_len, 1, n_land_tile, nsoilt,          &
                         sm_levels, nsmax, nice_use, n_wtrac_jls, n_evap_srce, &
                         river_row_length_dum, river_rows_dum, l_wtrac_jls,    &
                         wtrac_jls_data)
    call wtrac_jls_assoc(wtrac_jls, wtrac_jls_data)

    ! Note, jules_vars needs setting up after the change to pdims_s below so is
    ! not with the rest of these allocations.

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    !-----------------------------------------------------------------------
    ! surface forcing
    if ( iseasurfalg == iseasurfalg_specified_roughness ) then
      l_spec_z0 = .true.
      z0m_specified = z0m_specified_nml
      z0h_specified = z0h_specified_nml
      do i = 1, seg_len
        z0msea_2d(map_2d(1,i)) = z0m_specified(i,1)
      end do
    else
      l_spec_z0 = .false.
    end if

    ! Size this with stencil for use in Jules routines called
    pdims_s%i_start=1
    pdims_s%i_end=seg_len_halo

    ! Initialise those fields whose contents will not be fully set
    allocate(fracaero_t(land_field,ntiles))
    allocate(fracaero_s(land_field,ntiles))
    allocate(epot_surft(land_field,ntiles))
    do n = 1, ntiles
      do l = 1, land_field
        fracaero_t(l,n) = 0.0_r_um
        fracaero_s(l,n) = 0.0_r_um
        epot_surft(l,n) = 0.0_r_um
      end do
    end do

    call jules_vars_alloc(land_field,ntype,n_land_tile,rad_nband,nsoilt,       &
                          sm_levels, seg_len, 1, npft, bl_levels, pdims_s,     &
                          pdims, l_albedo_obs, cansnowtile, l_deposition,        &
                          jules_vars_data)
    call jules_vars_assoc(jules_vars,jules_vars_data)

    if (can_rad_mod == 6) then
      jules_vars%diff_frac = 0.4_r_um
    end if

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into Jules variables
    !-----------------------------------------------------------------------
    do i = 1, seg_len_halo
      k_blend_uv(i,1) = 1
    end do
    do i = 1, seg_len
      k_blend_tq(i,1) = 1
    end do

    l = 0
    do i = 1, seg_len
      if (flandg(i,1) > 0.0_r_um) then
        l = l+1
        ainfo%land_index(l) = i
      end if
    end do

    ! Sea-ice fraction
    do i = 1, seg_len
      do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = n - first_sea_ice_tile + 1
        ainfo%ice_fract_ij(i,1) = ainfo%ice_fract_ij(i,1) + &
             real(tile_fraction(tile_stencil(1,1,i)+n-1), r_um)
        ainfo%ice_fract_ncat_sicat(i, 1, i_sice) = &
             real(tile_fraction(tile_stencil(1,1,i)+n-1), r_um)
      end do

      ! Because Jules tests on flandg < 1, we need to ensure this is exactly
      ! 1 when no sea or sea-ice is present
      if ( tile_fraction(tile_stencil(1,1,i)+first_sea_tile-1) == 0.0_r_def .and. &
           ainfo%ice_fract_ij(i,1) == 0.0_r_um) then
        flandg(i,1) = 1.0_r_um
      end if

      ! Jules requires sea-ice fractions with respect to the sea area
      if (ainfo%ice_fract_ij(i, 1) > 0.0_r_um) then
        ainfo%ice_fract_ij(i, 1) = ainfo%ice_fract_ij(i, 1) / (1.0_r_um - flandg(i, 1))
        ainfo%ice_fract_ncat_sicat(i, 1, 1:n_sea_ice_tile) &
           = ainfo%ice_fract_ncat_sicat(i, 1, 1:n_sea_ice_tile) / (1.0_r_um - flandg(i, 1))
      end if
    end do

    do l = 1, land_field
      coast%fland(l) = flandg(ainfo%land_index(l),1)
      do n = 1, n_land_tile
        ! Jules requires fractions with respect to the land area
        ainfo%frac_surft(l, n) = real(tile_fraction(tile_stencil(1,1,ainfo%land_index(l))+n-1), r_um) &
             / coast%fland(l)
      end do
    end do

    if (buddy_sea == buddy_sea_on) then
      ! Set up map fields
      do i = 1, seg_len
        do m = 1, tile_stencil_size(i)
          coast%mapi(i,m) = tile_stencil(1,m,i)/n_surf_tile+1
        end do
        ! Deal with corners of cube sphere
        if (tile_stencil_size(i) == 8) then
          ! If the point doesn't exist, set it to the current central point
          ! since this is guaranteed not to be a sea point if it is a coastal
          ! point, hence will be excluded from the averaging
          coast%mapi(i,9) = tile_stencil(1,1,i)/n_surf_tile+1
        end if
      end do
      coast%mapj(1,1) = 1
    else
      do i = 1, seg_len_halo
        flandfac(i,1) = 1.0_r_def
        fseafac(i,1) = 1.0_r_def
      end do
    end if
    do i = 1, seg_len_halo
      jules_vars%u_1_p_ij(i,1) = u_in_w3(u_w3_stencil(1,1,i)+k_blend_uv(i,1)-1)
      jules_vars%v_1_p_ij(i,1) = v_in_w3(v_w3_stencil(1,1,i)+k_blend_uv(i,1)-1)
    end do

    ! Surface currents
    ! Coupled A-O models will have ocean surface current components on w3
    ! available after coupling exchanges, assuming that coupled models always
    ! couple to an ocean model AND that sea surface currents are always exchanged.
    ! If sea surface currents are obtained from elsewhere, e.g. ancillary files,
    ! then we need to ensure that sea_u_current and sea_v_current are populated
    ! from the appropriate source data.
    do i = 1, seg_len_halo
      jules_vars%u_0_p_ij(i,1) = sea_u_current(sea_u_w3_stencil(1,1,i))
      jules_vars%v_0_p_ij(i,1) = sea_v_current(sea_v_w3_stencil(1,1,i))
    end do

    ! Set type_pts and type_index
    call tilepts(land_field, ainfo%frac_surft, ainfo%surft_pts,                &
                 ainfo%surft_index, ainfo%l_lice_point, ainfo%l_lice_surft)

    ! combined sea and sea-ice index
    ssi_pts = seg_len
    do i = 1, seg_len
      if (flandg(i, 1) < 1.0_r_um) then
        ainfo%ssi_index(i) = i
      end if
      ainfo%fssi_ij(i,1) = 1.0_r_um - flandg(i, 1)
    end do

    ! individual sea and sea-ice indices
    ! first set defaults
    sea_pts = 0
    ! then calculate based on state
    do i = 1, seg_len
      if (ainfo%ssi_index(i) > 0) then
        if (ainfo%ice_fract_ij(i, 1) < 1.0_r_um) then
          sea_pts = sea_pts + 1
          ainfo%sea_index(sea_pts) = i
          ainfo%sea_frac(i) = 1.0_r_um - ainfo%ice_fract_ij(i,1)
        end if
      end if
    end do

    ! multi-category sea-ice index
    do i = 1, seg_len
      do n = 1, nice_use
        if (ainfo%ssi_index(i) > 0 .and. ainfo%ice_fract_ncat_sicat(i, 1, n) > 0.0_r_um) then
          ainfo%sice_pts_ncat(n) = ainfo%sice_pts_ncat(n) + 1
          ainfo%sice_index_ncat(ainfo%sice_pts_ncat(n), n) = i
          ainfo%sice_frac_ncat(i, n) = ainfo%ice_fract_ncat_sicat(i, 1, n)
        end if
      end do
    end do

    ! Land tile temperatures
    do i = 1, seg_len
      tstar_land(i,1) = 0.0_r_um
    end do
    do l = 1, land_field
      do n = 1, n_land_tile
        progs%tstar_surft(l, n) = real(tile_temperature(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        tstar_land(ainfo%land_index(l),1) = tstar_land(ainfo%land_index(l),1) &
             + ainfo%frac_surft(l, n) * progs%tstar_surft(l, n)
      end do
    end do

    ! Sea temperature
    ! Default to temperature over frozen sea as the initialisation
    ! that follows does not initialise sea points if they are fully
    ! frozen
    do i = 1, seg_len
      if (tile_fraction(tile_stencil(1,1,i)+first_sea_tile-1) > 0.0_r_def) then
        coast%tstar_sea_ij(i,1) = real(tile_temperature(map_tile(1,i)+first_sea_tile-1), r_um)
      else
        coast%tstar_sea_ij(i,1) = tfs
      end if
    end do

    ! Sea-ice temperatures
    do i = 1, seg_len
      tstar_sice(i,1) = 0.0_r_um
      if (ainfo%ice_fract_ij(i, 1) > 0.0_r_um) then
        do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = n - first_sea_ice_tile + 1
          coast%tstar_sice_sicat(i, 1, i_sice) = real(tile_temperature(map_tile(1,i)+n-1), r_um)
          tstar_sice(i,1) = tstar_sice(i,1) &
                   + ainfo%ice_fract_ncat_sicat(i,1,i_sice) * &
                   coast%tstar_sice_sicat(i,1,i_sice) / ainfo%ice_fract_ij(i,1)
        end do
      end if
    end do

    do i = 1, seg_len
      ! Sea & Sea-ice temperature
      coast%tstar_ssi_ij(i,1) = (1.0_r_um - ainfo%ice_fract_ij(i,1)) * &
           coast%tstar_sea_ij(i,1) + ainfo%ice_fract_ij(i,1) * tstar_sice(i,1)

      ! Grid-box mean surface temperature
      fluxes%tstar_ij(i,1) = flandg(i,1) * tstar_land(i,1) &
           + (1.0_r_um - flandg(i,1)) * coast%tstar_ssi_ij(i,1)
    end do

    ! Sea-ice conductivity and bulk temperature
    do i = 1, seg_len
      if (ainfo%ice_fract_ij(i, 1) > 0.0_r_um) then
        do n = 1, n_sea_ice_tile
          progs%k_sice_sicat(i, 1, n) = real(sea_ice_conductivity(map_sice(1,i)+n-1), r_um)
          ainfo%ti_cat_sicat(i, 1, n) = real(sea_ice_temperature(map_sice(1,i)+n-1), r_um)
          progs%ti_sicat(i,1,1) = progs%ti_sicat(i,1,1) &
                + ainfo%ice_fract_ncat_sicat(i,1,n) * ainfo%ti_cat_sicat(i,1,n) / ainfo%ice_fract_ij(i,1)
        end do
      end if
    end do

    do l = 1, land_field
      do n = 1, npft
        ! Leaf area index
        progs%lai_pft(l, n) = real(leaf_area_index(map_pft(1,ainfo%land_index(l))+n-1), r_um)
        ! Canopy height
        progs%canht_pft(l, n) = real(canopy_height(map_pft(1,ainfo%land_index(l))+n-1), r_um)
      end do
      ! Roughness length (z0_tile)
      psparms%z0m_soil_gb(l) = real(soil_roughness(map_2d(1,ainfo%land_index(l))), r_um)
    end do

    ! Urban ancillaries
    if ( l_urban2t ) then
      do l = 1, land_field
        urban_param%wrr_gb(l)   = real(urbwrr(map_2d(1,ainfo%land_index(l))), r_um)
        urban_param%hwr_gb(l)   = real(urbhwr(map_2d(1,ainfo%land_index(l))), r_um)
        urban_param%hgt_gb(l)   = real(urbhgt(map_2d(1,ainfo%land_index(l))), r_um)
        urban_param%ztm_gb(l)   = real(urbztm(map_2d(1,ainfo%land_index(l))), r_um)
        urban_param%disp_gb(l)  = real(urbdisp(map_2d(1,ainfo%land_index(l))), r_um)
      end do
    end if

    call sparm(land_field, n_land_tile, ainfo%surft_pts, ainfo%surft_index,   &
               ainfo%frac_surft, progs%canht_pft, progs%lai_pft,              &
               psparms%z0m_soil_gb, psparms%catch_snow_surft,                 &
               psparms%catch_surft, psparms%z0_surft, psparms%z0h_bare_surft, &
               urban_param%ztm_gb)


    do l = 1, land_field
      ! Snow-free soil albedo
      psparms%albsoil_soilt(l,1) = real(soil_albedo(map_2d(1,ainfo%land_index(l))), r_um)

      do m = 1, sm_levels
        ! Volumetric soil moisture at wilting point (smvcwt_soilt)
        psparms%smvcwt_soilt(l,1,m) = real(soil_moist_wilt(map_2d(1,ainfo%land_index(l))), r_um)
        ! Volumetric soil moisture at critical point (smvccl_soilt)
        psparms%smvccl_soilt(l,1,m) = real(soil_moist_crit(map_2d(1,ainfo%land_index(l))), r_um)
        ! Volumetric soil moisture at saturation (smvcst_soilt)
        psparms%smvcst_soilt(l,1,m) = real(soil_moist_sat(map_2d(1,ainfo%land_index(l))), r_um)
        ! Saturated soil water suction (sathh_soilt)
        psparms%sathh_soilt(l, 1, m) = real(soil_suction_sat(map_2d(1,ainfo%land_index(l))), r_um)
        ! Clapp and Hornberger b coefficient (bexp_soilt)
        psparms%bexp_soilt(l, 1, m) = real(clapp_horn_b(map_2d(1,ainfo%land_index(l))), r_um)
        ! Soil temperature (t_soil_soilt)
        progs%t_soil_soilt(l,1,m) = real(soil_temperature(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        ! Unfrozen soil moisture proportion (sthu_soilt)
        psparms%sthu_soilt(l,1,m) = real(unfrozen_soil_moisture(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        ! Frozen soil moisture proportion (sthf_soilt)
        psparms%sthf_soilt(l,1,m) = real(frozen_soil_moisture(map_soil(1,ainfo%land_index(l))+m-1), r_um)
      end do

      ! Soil thermal conductivity (hcon_soilt)
      psparms%hcon_soilt(l,1,:) = real(soil_thermal_cond(map_2d(1,ainfo%land_index(l))), r_um)

      ! Soil ancils dependant on smvcst_soilt (soil moisture saturation limit)
      if ( psparms%smvcst_soilt(l,1,1) > 0.0_r_um ) then
        ainfo%l_soil_point(l) = .true.
      end if
    end do

    ! Scaling factors needed for use in surface exchange code
    if (l_albedo_obs) then
      do l = 1, land_field
        if (coast%fland(l) > 0.0_r_um) then
          do n = 1, rad_nband
            do i = 1, n_land_tile
              i_tile=n_land_tile*(n-1) + i - 1
               jules_vars%albobs_scaling_surft(l,i,n) = &
                  albedo_obs_scaling(map_scal(1,ainfo%land_index(l))+i_tile)
            end do
          end do
        end if
      end do
    end if

    do i = 1, seg_len
      ! Cosine of the solar zenith angle
      psparms%cosz_ij(i,1) = real(cos_zen_angle(map_2d(1,i)), r_um)

      ! Downwelling LW radiation at surface
      forcing%lw_down_ij(i,1) = real(lw_down_surf(map_2d(1,i)), r_um)
    end do

    if (topography == topography_horizon) then
      ! Set skyview factor used internally by JULES
      do i = 1, seg_len
        sky(i,1) = real(skyview(map_2d(1,i)), r_um)
      end do
    end if

    ! Net SW radiation on tiles
    do l = 1, land_field
      do n = 1, n_land_tile
        ! Net SW radiation on tiles
        fluxes%sw_surft(l, n) = real(sw_down_surf(map_2d(1,ainfo%land_index(l))) - &
                                sw_up_tile(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        fluxes%emis_surft(l, n) = 1.0_r_um - &
             real(tile_lw_grey_albedo(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      end do
    end do

    ! Emissivity is set in lw_rad_tile_kernel
    do n = 1, n_land_tile
      fluxes%l_emis_surft_set(n) = .true.
    end do

    do i = 1, seg_len
      ! Net SW on open sea
      fluxes%sw_sea(i) = real(sw_down_surf(map_2d(1,i)) - &
                         sw_up_tile(map_tile(1,i)+first_sea_tile-1), r_um)

      ! The amount of diffuse visible light is needed for solar penetrating radiation
      sw_diffuse_blue_surf = sw_down_blue_surf(map_2d(1,i)) - sw_direct_blue_surf(map_2d(1,i))

      ! Calculate penetrating solar into sea ice as a combination of direct and diffuse
      ! visible light, multipled by the fraction penetrating
      if (l_sice_swpen .and. ocn_cpl_point(map_2d(1,i)) == 1_i_def) then
        do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = n - first_sea_ice_tile + 1
          sea_ice_pensolar(map_sice(1,i)+i_sice-1) = sw_direct_blue_surf(map_2d(1,i)) * &
               sea_ice_pensolar_frac_direct(map_sice(1,i)+i_sice-1) + &
               sw_diffuse_blue_surf *                                 &
               sea_ice_pensolar_frac_diffuse(map_sice(1,i)+i_sice-1)
        end do
      else
        do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = n - first_sea_ice_tile + 1
          sea_ice_pensolar(map_sice(1,i)+i_sice-1) = 0.0_r_def
        end do
      endif

      do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = n - first_sea_ice_tile + 1

        ! Net SW on sea-ice
        fluxes%sw_sicat(i, i_sice) = real(sw_down_surf(map_2d(1,i)) -  &
                                     sw_up_tile(map_tile(1,i)+n-1)  -  &
                                     sea_ice_pensolar(map_sice(1,i)+i_sice-1), r_um)
      end do

      ! photosynthetically active downwelling SW radiation
      photosynth_act_rad(i,1) = real(sw_down_blue_surf(map_2d(1,i)), r_um)
    end do

    do l = 1, land_field
      ! Ozone
      chemvars%o3_gb(l) = real(ozone(map_wth(1,ainfo%land_index(l))), r_um)
    end do

    ! Carbon dioxide
    co2_mmr = real(co2_mix_ratio_now, r_um)
    co2 = co2_mmr

    do l = 1, land_field
      ! Half of peak-to-trough height over root(2) of orography (ho2r2_orog_gb)
      jules_vars%ho2r2_orog_gb(l) = real(peak_to_trough_orog(map_2d(1,ainfo%land_index(l))), r_um)

      jules_vars%sil_orog_land_gb(l) = real(silhouette_area_orog(map_2d(1,ainfo%land_index(l))), r_um)

      ! Surface conductance (gs_gb)
      progs%gs_gb(l) = real(surface_conductance(map_2d(1,ainfo%land_index(l))), r_um)

      ! Canopy water on each tile (canopy_surft)
      do n = 1, n_land_tile
        progs%canopy_surft(l, n) = real(canopy_water(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      end do

      i_snow = 0
      do n = 1, n_land_tile
        progs%snow_surft(l,n) = real(tile_snow_mass(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        progs%nsnow_surft(l,n) = n_snow_layers(map_tile(1,ainfo%land_index(l))+n-1)
        progs%snowdepth_surft(l,n) = real(snow_depth(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        do j = 1, nsmax
          progs%ds_surft(l,n,j) = real(snow_layer_thickness(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          progs%sice_surft(l,n,j) = real(snow_layer_ice_mass(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          progs%sliq_surft(l,n,j) = real(snow_layer_liq_mass(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          progs%tsnow_surft(l,n,j) = real(snow_layer_temp(map_snow(1,ainfo%land_index(l))+i_snow), r_um)
          ! Counting from 0 so increment here
          i_snow = i_snow + 1
        end do
      end do

    end do

    ! Dust fields
    do i = 1, seg_len
      soil_clay(i, 1) = real(soil_clay_2d(map_2d(1,i)), r_um)
    end do
    do k = 1, ndiv
      do i = 1, seg_len
        dust_flux(i,1,k) = 0.0_r_um
      end do
    end do

    ! Level heights
    do i = 1, seg_len
      ainfo%z1_tq_ij(i,1) = height_wth(map_wth(1,i) + k_blend_tq(i,1))  &
                          - height_wth(map_wth(1,i))
      ainfo%z1_uv_ij(i,1) = height_w3(map_w3(1,i) + k_blend_uv(i,1)-1)  &
                          - height_wth(map_wth(1,i))
      z1_uv_top(i,1) = height_wth(map_wth(1,i) + 1) - height_wth(map_wth(1,i))
      z1_tq_top(i,1) = height_w3(map_w3(1,i) + 1) - height_wth(map_wth(1,i))
    end do

    l=0
    do i = 1, seg_len
      ! Land height
      if (land_sea_mask(i,1)) then
        l=l+1
        if ( (l_ctile .and. coast%fland(l) > 0.0_r_um .and. &
                            coast%fland(l) < 1.0_r_um) &
                            .or. any(l_elev_absolute_height) ) then
          jules_vars%z_land_ij(i,1) = (height_wth(map_wth(1,i))+planet_radius)-planet_radius
          if (jules_vars%z_land_ij(i,1) <  0.0_r_um) &
              jules_vars%z_land_ij(i,1) = 0.0_r_um
        end if
      end if
    end do

    do i = 1, seg_len
      ! thermodynamic variables
      temperature(i,1,1) = theta_in_wth(map_wth(1,i)+k_blend_tq(i,1)) * &
                         exner_in_wth(map_wth(1,i)+k_blend_tq(i,1))
      q(i,1,1) = m_v_n(map_wth(1,i)+k_blend_tq(i,1))
      qcl(i,1,1) = m_cl_n(map_wth(1,i)+k_blend_tq(i,1))
      if (l_noice_in_turb) then
        qcf(i,1,1) = 0.0_r_um
        bulk_cloud_fraction(i,1,1) = cf_liquid(map_wth(1,i)+k_blend_tq(i,1))
      else
        qcf(i,1,1) = m_cf_n(map_wth(1,i)+k_blend_tq(i,1))
        bulk_cloud_fraction(i,1,1) = cf_bulk(map_wth(1,i)+k_blend_tq(i,1))
      end if
      forcing%qw_1_ij(i,1) = q(i,1,1) + qcl(i,1,1) + qcf(i,1,1)
      forcing%tl_1_ij(i,1) = temperature(i,1,1) - lcrcp*qcl(i,1,1) - lsrcp*qcf(i,1,1)

      ! pressure
      p_theta_levels(i,1,1) = p_zero*(exner_in_wth(map_wth(1,i)+k_blend_tq(i,1)))**(1.0_r_def/kappa)
      forcing%pstar_ij(i,1) = p_zero*(exner_in_wth(map_wth(1,i) + 0))**(1.0_r_def/kappa)
    end do

    !-----------------------------------------------------------------------
    ! Things saved from one timestep to the next
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      ! previous BL height
      jules_vars%zh(i,1) = zh_2d(map_2d(1,i))
      ! surface roughness
      progs%z0msea_ij(i,1) = z0msea_2d(map_2d(1,i))
      ! downdraft at cloud base
      ddmfx(i,1) = dd_mf_cb(map_2d(1,i))
    end do
    ! needed to ensure neutral diagnostics can be calculated
    sf_diag%suv10m_n  = .true.
    ! needed to ensure z0h_eff is saved if wanted
    sf_diag%l_z0h_eff_gb = .not. associated(z0h_eff, empty_real_data)
    sf_diag%l_z0m_gb = .true.
    call alloc_sf_expl(sf_diag, .true., land_field)

    !-----------------------------------------------------------------------
    ! External science code called
    !-----------------------------------------------------------------------

    call buoy_tq (                                                             &
       ! IN dimensions/logicals
       1,                                                                      &
       ! IN fields
       p_theta_levels,temperature,q,qcf,qcl,bulk_cloud_fraction,               &
       ! OUT fields
       bt,bq,bt_cld,bq_cld,bt_blend,bq_blend,a_qs,a_dqsdt,dqsdt                &
       )

    allocate(hcons_soilt(land_field))
    allocate(emis_soil(land_field))
    allocate(dtstar_surft(land_field,ntiles))
    allocate(alpha1(land_field,ntiles))
    allocate(ashtf_prime_surft(land_field,ntiles))
    allocate(chr1p5m(land_field,ntiles))
    allocate(flake(land_field,ntiles))
    allocate(resft(land_field,ntiles))
    allocate(resfs(land_field,ntiles))
    allocate(rhokh_surft(land_field,ntiles))
    allocate(canhc_surft(land_field,ntiles))
    allocate(wt_ext_surft(land_field,sm_levels,ntiles))
    allocate(resp_s_tot_soilt(land_field))
    call surf_couple_explicit(                                                 &
      !Arguments used by JULES-standalone
      !Misc INTENT(IN) DONE
      bq_blend, bt_blend, photosynth_act_rad, day_of_year,                     &
      second_of_day,                                                           &
      !INOUT Diagnostics, in sf_diags_mod
      sf_diag,                                                                 &
      !Fluxes INTENT(OUT) DONE
      fqw,ftl,                                                                 &
      !Misc INTENT(OUT)
      !rhokms needed for message passing
      radnet_sice, rhokm, rhokm_land, rhokm_ssi,                               &
      !Out of explicit and into implicit only INTENT(OUT)
      !cdr10m needed for message passing
      cdr10m, alpha1, alpha1_sea, alpha1_sice, ashtf_prime, ashtf_prime_sea,   &
      ashtf_prime_surft, epot_surft,                                           &
      !rhokh needed in BL
      fracaero_t, fracaero_s, resfs, resft,                                    &
      rhokh, rhokh_surft, rhokh_sice, rhokh_sea,                               &
      dtstar_surft, dtstar_sea, dtstar_sice, z0hssi, z0mssi, chr1p5m,          &
      chr1p5m_sice, canhc_surft, wt_ext_surft, flake,                          &
      !Out of explicit and into extra only INTENT(OUT)
      hcons_soilt,                                                             &
      !Out of explicit and into implicit and extra INTENT(OUT)
      ainfo%frac_surft,                                                        &
      !Additional arguments for the BL-----------------------------------------
      !JULES prognostics module
      !IN
      ainfo%ti_cat_sicat,                                                      &
      !JULES ancil_info module
      !IN
      land_field, ssi_pts, sea_pts, ntiles,                                    &
      !OUT
      ainfo%surft_pts,                                                         &
      ! IN input data from the wave model
      charnock_w,                                                              &
      !JULES coastal module
      !4 IN, vshr_ both OUT
      flandg,                                                                  &
      !IN
      co2,                                                                     &
      !JULES trifctl module
      !IN
      asteps_since_triffid,                                                    &
      !JULES p_s_parms module
      !IN
      soil_clay,                                                               &
      !JULES switches module  **squish**
      !IN SCM related
      l_spec_z0,                                                               &
      !Not in a JULES module
      !IN
      1, 1, z1_uv_top, z1_tq_top, sky, ddmfx,                                  &
      !3 IN, 1 OUT requiring STASH flag
      l_aero_classic, z0m_specified, z0h_specified,                            &
      !OUT not requiring STASH flag
      recip_l_mo_sea, rib_gb,                                                  &
      !OUT 2 message passing, 1 soil moisture nudging, rest of BL
      flandfac, fseafac, fb_surf, ustargbm, t1_sd, q1_sd, rhostar,             &
      !OUT
      vshr, resp_s_tot_soilt, emis_soil,                                       &
      !TYPES containing field data (IN OUT)
      crop_vars,psparms,ainfo,trif_vars,aerotype,urban_param,progs,trifctltype,&
      coast, jules_vars, fluxes, lake_vars, forcing, chemvars, wtrac_jls,      &
      progs_cbl_vars, work_cbl                                                 &
      )

    if (flux_bc_opt > interactive_fluxes) then
      ! Surface fluxes calculated in SFEXPL are substituted with
      ! forcing values.  NOTE: Surface calculation also made
      ! explicit (time weight set to zero).
      if (flux_bc_opt == specified_fluxes_cd) then
        do i = 1, seg_len
          ! In the future we will set rhokm using input ustar_in.
          ! For now, we set using a fixed value.
          uStarGBM(i,1) = ustar_fixed_value
          rhokm(i,1) = rhostar(i,1)*uStarGBM(i,1)*uStarGBM(i,1)/vshr(i,1)
          if (flandg(i,1) > 0.0_r_um) then
            rhokm_land(i,1) = rhostar(i,1)*uStarGBM(i,1)*uStarGBM(i,1)/        &
                              coast%vshr_land_ij(i,1)
          end if
          if (flandg(i,1) < 1.0_r_um) then
            rhokm_ssi(i,1)  = rhostar(i,1)*uStarGBM(i,1)*uStarGBM(i,1)/        &
                              coast%vshr_ssi_ij(i,1)
          end if
        end do
      end if

      do i = 1, seg_len
        !..Converts Fluxes from W/m^2 to rho*K/s
        rholem = rhostar(i,1)

        !..If comparing against LES with rho ne rhostar then match
        ! w'theta' rather than rho*wtheta
        ! RHOLEM = 1.0

        fqw(i,1)   = (rhostar(i,1)*flux_e)/(lc*rholem)
        ftl(i,1)   = (rhostar(i,1)*flux_h)/(cp*rholem)

        fb_surf(i,1) = g * ( bt_blend(i,1,1)*ftl(i,1) +                        &
                             bq_blend(i,1,1)*fqw(i,1) ) /rhostar(i,1)
        recip_l_mo_sea(i,1) = -vkman * fb_surf(i,1)                           &
                              / ( ustargbm(i,1)*ustargbm(i,1)*ustargbm(i,1) )
        ! Zeroing w_m and tv1_sd required for Psyclone transmute
        w_m = 0.0_r_um
        tv1_sd = 0.0_r_um
        if ( fb_surf(i,1)  >   0.0_r_um) then
          w_m  = ( 0.25_r_um*jules_vars%zh(i,1)*fb_surf(i,1) +                &
                  ustargbm(i,1)*ustargbm(i,1)*ustargbm(i,1) ) ** one_third

          t1_sd(i,1) = 1.93_r_um * ftl(i,1) / (rhostar(i,1) * w_m)
          q1_sd(i,1) = 1.93_r_um * fqw(i,1) / (rhostar(i,1) * w_m)
          tv1_sd     = temperature(i,1,1) * ( bt_blend(i,1,1)*t1_sd(i,1) +     &
                                            bq_blend(i,1,1)*q1_sd(i,1) )
          t1_sd(i,1) = max ( 0.0_r_um , t1_sd(i,1) )
          q1_sd(i,1) = max ( 0.0_r_um , q1_sd(i,1) )
          if (tv1_sd  <=  0.0_r_um) then
            t1_sd(i,1) = 0.0_r_um
            q1_sd(i,1) = 0.0_r_um
          end if
        else
          t1_sd(i,1) = 0.0_r_um
          q1_sd(i,1) = 0.0_r_um
        end if
      end do
      do l = 1, land_field
        do n = 1, ntiles
          fluxes%fqw_surft(l,n) = fqw(ainfo%land_index(l),1)
          fluxes%ftl_surft(l,n) = ftl(ainfo%land_index(l),1)
        end do
      end do

      if ( flux_bc_opt == specified_fluxes_only ) then
        ! recalculate the surface temperature here to be consistent with
        ! fixed fluxes, since this isn't done in conv_surf_flux or specified
        allocate(qs_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
        do i = 1, seg_len
          ! start with simple extrapolation from level 1
          fluxes%tstar_ij(i,1) = temperature(i,1,1) + grcp * ainfo%z1_tq_ij(i,1)

          call qsat_mix(qs_star,fluxes%tstar_ij,forcing%pstar_ij,pdims%i_len,pdims%j_len)

          dqsdt_star = repsilon * lc * qs_star(i,1) /                          &
                       ( r * fluxes%tstar_ij(i,1) * fluxes%tstar_ij(i,1) )

          theta1 = temperature(i,1,1) * (p_zero/p_theta_levels(i,1,1))**kappa

          wthvbar = theta1 *                                                   &
                    (1.0_r_um+c_virtual*q(i,1,1)-qcl(i,1,1)-qcf(i,1,1)) *      &
                    fb_surf(i,1) / g

          ch = rhokh(i,1) / ( vshr(i,1) * rhostar(i,1) )

          ! Now more complicated formula based on fluxes
          fluxes%tstar_ij(i,1) = ( theta1 + wthvbar/(ch*max(0.1_r_um,vshr(i,1))) -  &
               c_virtual * theta1 *                                            &
               (qs_star(i,1)-q(i,1,1)-dqsdt_star*fluxes%tstar_ij(i,1)) )       &
               / ( (p_zero/forcing%pstar_ij(i,1))**kappa +                     &
               c_virtual * theta1 * dqsdt_star )

          tstar_land(i,1) = fluxes%tstar_ij(i,1)
          coast%tstar_sea_ij(i,1)  = fluxes%tstar_ij(i,1)
          coast%tstar_sice_sicat(i,1,:) = fluxes%tstar_ij(i,1)
          coast%tstar_ssi_ij(i,1)  = fluxes%tstar_ij(i,1)

        end do
        deallocate(qs_star)

      end if
    end if ! flux_bc_opt

    !-----------------------------------------------------------------------
    ! Mineral dust production
    !-----------------------------------------------------------------------
    if (l_dust_flux_only) then
      allocate(dust_flux_surft(land_field,ntiles,ndiv))
      allocate(u_s_t_tile(land_field,ntiles,ndivh))
      allocate(u_s_t_dry_tile(land_field,ntiles,ndivh))
      allocate(soil_layer_moisture(land_field,sm_levels))
      allocate(rhostar_land(land_field))
      allocate(sand_land(land_field))
      allocate(clay_land(land_field))
      allocate(mrel_land(land_field,ndivl))
      allocate(dust_emsc_land(land_field))
      !put fields into land arrays
      do l = 1, land_field
        rhostar_land(l) = rhostar(ainfo%land_index(l),1)
        sand_land(l)   = real(soil_sand_2d(map_2d(1,ainfo%land_index(l))), r_um)
        clay_land(l)   = real(soil_clay_2d(map_2d(1,ainfo%land_index(l))), r_um)
        mrel_land(l,1) = real(dust_div_mrel(map_dust(1,ainfo%land_index(l)) + 0), r_um)
        mrel_land(l,2) = real(dust_div_mrel(map_dust(1,ainfo%land_index(l)) + 1), r_um)
        mrel_land(l,3) = real(dust_div_mrel(map_dust(1,ainfo%land_index(l)) + 2), r_um)
        mrel_land(l,4) = real(dust_div_mrel(map_dust(1,ainfo%land_index(l)) + 3), r_um)
        mrel_land(l,5) = real(dust_div_mrel(map_dust(1,ainfo%land_index(l)) + 4), r_um)
        mrel_land(l,6) = real(dust_div_mrel(map_dust(1,ainfo%land_index(l)) + 5), r_um)
        dust_emsc_land(l) = 1.0_r_um
        ! Soil moisture content (kg m-2, soil_layer_moisture)
        do m = 1, sm_levels
          soil_layer_moisture(l,m) = real(soil_moisture(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        end do
      end do !land_field

      call dust_srce(                                                          &
           ! IN arguments
             land_field,ntiles,ainfo%surft_pts,ainfo%surft_index,coast%fland,  &
             progs%tstar_surft,rhostar_land,soil_layer_moisture,               &
             progs%snow_surft, aerotype%u_s_std_surft,mrel_land,clay_land,     &
             sand_land, jules_vars%ho2r2_orog_gb,dust_emsc_land,               &
             ! OUT arguments
             dust_flux_surft,u_s_t_tile,u_s_t_dry_tile                         &
             )

      allocate(dust_emiss_frac(land_field,ntiles))
      ! Get the fraction within each tile which is bare soil, for the purpose
      ! of dust emission:
      call dust_calc_emiss_frac(                                               &
         land_field,ntiles,ainfo%surft_pts,ainfo%surft_index,ainfo%frac_surft, &
         progs%lai_pft, dust_emiss_frac                                        &
         )

      ! Produce a total dust flux over all tiles, by looping through tiles and
      ! multiplying the flux on that tile by the dust_emiss_frac, and summing.
      ! Note that in the following loop each land point i,j may have
      ! multiple tiles. Therefore, a COLLAPSE(2) clause will introduce a data
      ! race on dust_flux.
      do idiv = 1, ndiv
        do m = 1, ntiles
          do n = 1, ainfo%surft_pts(m)
            l = ainfo%surft_index(n,m)
            j = (ainfo%land_index(l)-1)/pdims%i_end + 1
            i = ainfo%land_index(l) - (j-1)*pdims%i_end
            dust_flux(i,j,idiv) = dust_flux(i,j,idiv) +                        &
                 dust_flux_surft(l,m,idiv)*dust_emiss_frac(l,m)
          end do !surft_pts
        end do !ntiles
      end do !ndiv

      deallocate(dust_emiss_frac)
      deallocate(mrel_land)
      deallocate(clay_land)
      deallocate(sand_land)
      deallocate(rhostar_land)
      deallocate(soil_layer_moisture)
      deallocate(u_s_t_dry_tile)
      deallocate(u_s_t_tile)
      deallocate(dust_flux_surft)

    end if !l_dust_flux

    do i = 1, seg_len
      ! variables passed to explicit BL
      rhostar_2d(map_2d(1,i)) = rhostar(i,1)
      recip_l_mo_sea_2d(map_2d(1,i)) = recip_l_mo_sea(i,1)
      t1_sd_2d(map_2d(1,i)) = t1_sd(i,1)
      q1_sd_2d(map_2d(1,i)) = q1_sd(i,1)

      ! 2D variables that need interpolating to cell faces
      surf_interp(map_surf(1,i)+0) = flandg(i,1)
      surf_interp(map_surf(1,i)+1) = rhokm_land(i,1)
      surf_interp(map_surf(1,i)+2) = rhokm_ssi(i,1)
      surf_interp(map_surf(1,i)+3) = flandfac(i,1)
      surf_interp(map_surf(1,i)+4) = fseafac(i,1)
      surf_interp(map_surf(1,i)+6) = fb_surf(i,1)
      surf_interp(map_surf(1,i)+7) = real(k_blend_uv(i,1)-1, r_def)
      surf_interp(map_surf(1,i)+8) = cdr10m(i,1)
      surf_interp(map_surf(1,i)+9) = real(sf_diag%cd10m_n(i,1), r_def)
      surf_interp(map_surf(1,i)+10)= real(sf_diag%cdr10m_n(i,1), r_def)

      ! surface forcing of atmosphere
      rhokm_bl(map_wth(1,i)) = rhokm(i,1)
      rhokh_bl(map_w3(1,i)) = rhokh(i,1)
      moist_flux_bl(map_w3(1,i)) = fqw(i,1)
      heat_flux_bl(map_w3(1,i)) = ftl(i,1)
      gradrinr(map_wth(1,i)) = rib_gb(i,1)
    end do

    ! surface prognostics and things needed elsewhere
    do n = 1, n_land_tile
      do l = 1, land_field
        alpha1_tile(map_tile(1,ainfo%land_index(l))+n-1) = alpha1(l,n)
        fracaero_t_tile(map_tile(1,ainfo%land_index(l))+n-1) = fracaero_t(l,n)
        fracaero_s_tile(map_tile(1,ainfo%land_index(l))+n-1) = fracaero_s(l,n)
        z0h_tile(map_tile(1,ainfo%land_index(l))+n-1) = fluxes%z0h_surft(l,n)
        z0m_tile(map_tile(1,ainfo%land_index(l))+n-1) = fluxes%z0m_surft(l,n)
        chr1p5m_tile(map_tile(1,ainfo%land_index(l))+n-1) = chr1p5m(l,n)
        resfs_tile(map_tile(1,ainfo%land_index(l))+n-1) = resfs(l,n)
        gc_tile(map_tile(1,ainfo%land_index(l))+n-1) = progs%gc_surft(l,n)
        canhc_tile(map_tile(1,ainfo%land_index(l))+n-1) = canhc_surft(l,n)
      end do
      ! Fields only calculated on tiles which exist
      do m = 1, ainfo%surft_pts(n)
        l = ainfo%surft_index(m,n)
        ashtf_prime_tile(map_tile(1,ainfo%land_index(l))+n-1) = ashtf_prime_surft(l,n)
        dtstar_tile(map_tile(1,ainfo%land_index(l))+n-1) = dtstar_surft(l,n)
        rhokh_tile(map_tile(1,ainfo%land_index(l))+n-1) = rhokh_surft(l,n)
      end do
    end do

    do n = 1, n_land_tile
      do m = 1, sm_levels
        i_tile=sm_levels*(n-1)+ m - 1 
        do k = 1, ainfo%surft_pts(n)
          l = ainfo%surft_index(k,n)
          tile_water_extract(map_smtile(1,ainfo%land_index(l))+i_tile) = wt_ext_surft(l,m,n)
        end do
      end do
    end do

    do i = 1, seg_len
      alpha1_tile(map_tile(1,i)+first_sea_tile-1) = alpha1_sea(i,1)
      z0m_tile(map_tile(1,i)+first_sea_tile-1) = progs%z0msea_ij(i,1)
    end do
    ! These are only calculated on sea points
    do l = 1, sea_pts
      i = ainfo%sea_index(l)
      ashtf_prime_tile(map_tile(1,i)+first_sea_tile-1) = ashtf_prime_sea(i,1)
      dtstar_tile(map_tile(1,i)+first_sea_tile-1) = dtstar_sea(i,1)
    end do
    ! This is only used with multi-category sea-ice
    if (n_sea_ice_tile > 1) then
      do i = 1, seg_len
        rhokh_tile(map_tile(1,i)+first_sea_tile-1) = rhokh_sea(i,1)
      end do
    end if

    do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = n - first_sea_ice_tile + 1
      do i = 1, seg_len
        alpha1_tile(map_tile(1,i)+n-1) = alpha1_sice(i,1,i_sice)
        rhokh_tile(map_tile(1,i)+n-1) = rhokh_sice(i,1,i_sice)
      end do
      ! Fields are only written on sea-ice points
      do l = 1, ainfo%sice_pts_ncat(i_sice)
        i = ainfo%sice_index_ncat(l,i_sice)
        ashtf_prime_tile(map_tile(1,i)+n-1) = ashtf_prime(i,1,i_sice)
        dtstar_tile(map_tile(1,i)+n-1) = dtstar_sice(i,1,i_sice)
      end do
    end do

    do i = 1, seg_len
      z0h_tile(map_tile(1,i)+first_sea_ice_tile-1) = z0hssi(i,1)
      z0m_tile(map_tile(1,i)+first_sea_ice_tile-1) = z0mssi(i,1)
      chr1p5m_tile(map_tile(1,i)+first_sea_ice_tile-1) = chr1p5m_sice(i,1)

      blend_height_tq(map_2d(1,i)) = k_blend_tq(i,1)
      z0m_eff(map_2d(1,i)) = jules_vars%z0m_eff_ij(i,1)
      ustar(map_2d(1,i)) = ustargbm(i,1)
    end do

    do l = 1, land_field
      soil_moist_avail(map_2d(1,ainfo%land_index(l))) = progs%smc_soilt(l,1)
    end do

    do n = 1, npft
      do l = 1, land_field
        ! Unloading rate of snow from plant functional types
        snow_unload_rate(map_pft(1,ainfo%land_index(l))+n-1) = real(jules_vars%unload_backgrnd_pft(l, n), r_def)
      end do
    end do

    do n = 1, n_land_tile
      do l = 1, land_field
        ! Land tile temperatures
        tile_temperature(map_tile(1,ainfo%land_index(l))+n-1) = real(progs%tstar_surft(l,n), r_def)
        ! sensible heat flux
        tile_heat_flux(map_tile(1,ainfo%land_index(l))+n-1) = real(fluxes%ftl_surft(l,n), r_def)
        ! moisture flux
        tile_moisture_flux(map_tile(1,ainfo%land_index(l))+n-1) = real(fluxes%fqw_surft(l,n), r_def)
      end do
    end do

    do i = 1, seg_len
      ! Sea temperature
      tile_temperature(map_tile(1,i)+first_sea_tile-1) = real(coast%tstar_sea_ij(i,1), r_def)
    end do

    do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = n - first_sea_ice_tile + 1
      do i = 1, seg_len
        ! sea-ice temperature
        tile_temperature(map_tile(1,i)+n-1) = real(coast%tstar_sice_sicat(i,1,i_sice), r_def)
        ! sea-ice heat flux
        tile_heat_flux(map_tile(1,i)+n-1) = real(fluxes%ftl_sicat(i,1,i_sice), r_def)
        ! sea-ice moisture flux
        tile_moisture_flux(map_tile(1,i)+n-1) = real(fluxes%fqw_sicat(i,1,i_sice), r_def)
      end do
    end do

    do i = 1, seg_len
      ! Cell roughness length for momentum
      z0m_2d(map_2d(1,i)) = real(sf_diag%z0m_gb(i,1), r_def)
      z0msea_2d(map_2d(1,i)) = progs%z0msea_ij(i,1)
    end do

    ! Dust fluxes
    do n = 1, ndiv
      do i = 1, seg_len
        dust_div_flux(map_dust(1,i) + n-1) = real(dust_flux(i,1,n), r_def)
      end do
    end do

    ! write diagnostics
    do l = 1, land_field
      net_prim_prod(map_2d(1,ainfo%land_index(l))) = real(trifctltype%npp_gb(l), r_def)
      surface_conductance(map_2d(1,ainfo%land_index(l))) = real(progs%gs_gb(l), r_def)
      thermal_cond_wet_soil(map_2d(1,ainfo%land_index(l))) = hcons_soilt(l)
    end do

    if (.not. associated(soil_respiration, empty_real_data) ) then
      if (dim_cs1 == 4) then
        do l = 1, land_field
          soil_respiration(map_2d(1,ainfo%land_index(l))) = resp_s_tot_soilt(l)
        end do
      else
        do l = 1, land_field
          soil_respiration(map_2d(1,ainfo%land_index(l))) = trifctltype%resp_s_soilt(l,1,1,1)
        end do
      end if
    end if

    if (.not. associated(gross_prim_prod, empty_real_data) ) then
      do i = 1, seg_len
        gross_prim_prod(map_2d(1,i)) = 0.0_r_def
      end do
      do l = 1, land_field
        gross_prim_prod(map_2d(1,ainfo%land_index(l))) = real(trifctltype%gpp_gb(l), r_def)
      end do
    end if

    if (.not. associated(z0h_eff, empty_real_data) ) then
      do i = 1, seg_len
        z0h_eff(map_2d(1,i)) = sf_diag%z0h_eff_gb(i,1)
      end do
    end if

    ! deallocate diagnostics
    call dealloc_sf_expl(sf_diag)

    ! Set back to SCM for use elsewhere
    pdims_s%i_start=1
    pdims_s%i_end=seg_len

    deallocate(resp_s_tot_soilt)
    deallocate(fracaero_t)
    deallocate(fracaero_s)
    deallocate(epot_surft)
    deallocate(hcons_soilt)
    deallocate(emis_soil)
    deallocate(dtstar_surft)
    deallocate(alpha1)
    deallocate(ashtf_prime_surft)
    deallocate(chr1p5m)
    deallocate(flake)
    deallocate(resft)
    deallocate(resfs)
    deallocate(rhokh_surft)
    deallocate(canhc_surft)
    deallocate(wt_ext_surft)

    call ancil_info_nullify(ainfo)
    call ancil_info_dealloc(ainfo_data)

    call forcing_nullify(forcing)
    call forcing_dealloc(forcing_data)

    call crop_vars_nullify(crop_vars)
    call crop_vars_dealloc(crop_vars_data)

    call lake_nullify(lake_vars)
    call lake_dealloc(lake_data)

    call coastal_nullify(coast)
    call coastal_dealloc(coastal_data)

    call trifctl_nullify(trifctltype)
    call trifctl_dealloc(trifctl_data)

    call urban_param_nullify(urban_param)
    call urban_param_dealloc(urban_param_data)

    call aero_nullify(aerotype)
    call aero_dealloc(aero_data)

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

    call chemvars_nullify(chemvars)
    call chemvars_dealloc(chemvars_data)

    call wtrac_jls_nullify(wtrac_jls)
    call wtrac_jls_dealloc(wtrac_jls_data)

  end subroutine jules_exp_code

end module jules_exp_kernel_mod
