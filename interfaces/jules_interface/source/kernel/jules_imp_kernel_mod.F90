!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the implicit UM boundary layer scheme.
!>
module jules_imp_kernel_mod

  use argument_mod,           only : arg_type,                  &
                                     GH_FIELD, GH_SCALAR,       &
                                     GH_INTEGER, GH_REAL,       &
                                     GH_READ, GH_WRITE, GH_INC, &
                                     GH_READWRITE, DOMAIN,      &
                                     ANY_DISCONTINUOUS_SPACE_1, &
                                     ANY_DISCONTINUOUS_SPACE_2, &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     ANY_DISCONTINUOUS_SPACE_4, &
                                     ANY_DISCONTINUOUS_SPACE_5, &
                                     ANY_DISCONTINUOUS_SPACE_6, &
                                     ANY_DISCONTINUOUS_SPACE_7
  use radiation_config_mod,      only : topography, topography_horizon
  use constants_mod,             only : i_def, i_um, r_def, r_um
  use empty_data_mod,            only : empty_real_data
  use fs_continuity_mod,         only : W3, Wtheta
  use kernel_mod,                only : kernel_type
  use timestepping_config_mod,   only : outer_iterations
  use water_constants_mod,       only : tfs, lc, lf

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: jules_imp_kernel_type
    private
    type(arg_type) :: meta_args(85) = (/                                          &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! outer
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! loop
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! wetrho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_fraction
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_thickness
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_pensolar
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_pensolar_frac_direct
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_pensolar_frac_diffuse
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! screen_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! time_since_transition
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! latitude
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_snow_mass
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! n_snow_layers
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! snow_depth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canopy_water
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! soil_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_heat_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_moisture_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! sw_up_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lw_down_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_blue_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_direct_blue_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! skyview
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_lw_grey_albedo
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! snowice_sublimation (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_heat_flux (W m-2)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! canopy_evap (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! water_extraction (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! snowice_melt (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! rh_crit_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rhokh_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! moist_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! heat_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! alpha1_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! ashtf_prime_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! dtstar_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! fracaero_t_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! fracaero_s_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! z0h_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! z0m_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! rhokh_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! chr1p5m_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! resfs_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canhc_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! tile_water_extract
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ustar
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lake_evap
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_avail
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_7),&! bl_type_ind
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTheta),                   &! qw_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTheta),                   &! tl_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! dqw1_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! dtl1_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ct_ctq1_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_ht_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! t1p5m_surft
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! q1p5m_surft
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! t1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! q1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qcl1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! rh1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! t1p5m_ssi
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! q1p5m_ssi
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qcl1p5m_ssi
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! rh1p5m_ssi
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! t1p5m_land
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! q1p5m_land
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qcl1p5m_land
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! rh1p5m_land
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! latent_heat
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! snomlt_surf_htf
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_evap
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_surf_ht_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_sw_net
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_radnet
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_lw_up
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_lw_down
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1) &! ocn_cpl_point
         /)

    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: jules_imp_code
  end type

  public :: jules_imp_code

contains

  !> @brief Interface to the implicit UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !>          NB This version uses winds in w3 space (i.e. A-grid)
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     seg_len              Number of horizontal points
  !> @param[in]     outer                Outer loop counter of dynamics SI scheme
  !> @param[in]     loop                 Loop counter of BL predictor-corrector scheme
  !> @param[in]     wetrho_in_w3         Wet density field in density space
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in]     tile_fraction        Surface tile fractions
  !> @param[in]     sea_ice_thickness    Depth of sea-ice (m)
  !> @param[in,out] sea_ice_temperature  Bulk temperature of sea-ice (K)
  !> @param[in,out] sea_ice_pensolar     Sea ice penetrating solar radiation (W m-2)
  !> @param[in] sea_ice_pensolar_frac_direct   Fraction of penetrating solar (direct visible) into sea ice
  !> @param[in] sea_ice_pensolar_frac_diffuse  Fraction of penetrating solar (diffuse visible) into sea ice
  !> @param[in,out] tile_temperature     Surface tile temperatures
  !> @param[in,out] screen_temperature   Tiled screen level liquid temperature
  !> @param[in,out] time_since_transition Time since decoupled screen transition
  !> @param[in]     latitude             Latitude of cell centre
  !> @param[in]     tile_snow_mass       Snow mass on tiles (kg/m2)
  !> @param[in]     n_snow_layers        Number of snow layers on tiles
  !> @param[in]     snow_depth           Snow depth on tiles
  !> @param[in]     canopy_water         Canopy water on each tile
  !> @param[in]     soil_temperature     Soil temperature
  !> @param[in,out] tile_heat_flux       Surface heat flux
  !> @param[in,out] tile_moisture_flux   Surface moisture flux
  !> @param[in]     sw_up_tile           Upwelling SW radiation on surface tiles
  !> @param[in]     sw_down_surf         Downwelling SW radiation at surface
  !> @param[in]     lw_down_surf         Downwelling LW radiation at surface
  !> @param[in]     sw_down_blue_surf    Photosynthetically active SW down
  !> @param[in]     sw_direct_blue_surf  Downwelling direct visible light
  !> @param[in]     skyview              Skyview / area enhancement factor
  !> @param[in]     tile_lw_grey_albedo  Surface tile longwave grey albedo
  !> @param[in,out] snowice_sublimation  Sublimation of snow and ice
  !> @param[in,out] surf_heat_flux       Surface heat flux
  !> @param[in,out] canopy_evap          Canopy evaporation from land tiles
  !> @param[in,out] water_extraction     Extraction of water from each soil layer
  !> @param[in,out] snowice_melt         Surface, canopy and sea ice, snow and ice melt rate
  !> @param[in]     m_cf                 Cloud frozen mixing ratio after advection
  !> @param[in]     rh_crit_wth          Critical relative humidity
  !> @param[in]     rhokh_bl             Heat eddy diffusivity on BL levels
  !> @param[in,out] moist_flux_bl        Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl         Vertical heat flux on BL levels
  !> @param[in]     dtrdz_tq_bl          dt/(rho*r*r*dz) in wth
  !> @param[in]     alpha1_tile          dqsat/dT in surface layer on tiles
  !> @param[in]     ashtf_prime_tile     Heat flux coefficient on tiles
  !> @param[in]     dtstar_tile          Change in surface temperature on tiles
  !> @param[in]     fracaero_t_tile      Fraction of moisture flux with only aerodynamic resistance
  !> @param[in]     fracaero_s_tile      Fraction of moisture flux with only aerodynamic resistance
  !> @param[in]     z0h_tile             Heat roughness length on tiles
  !> @param[in]     z0m_tile             Momentum roughness length on tiles
  !> @param[in]     rhokh_tile           Surface heat diffusivity on tiles
  !> @param[in]     chr1p5m_tile         1.5m transfer coefficients on tiles
  !> @param[in]     resfs_tile           Combined aerodynamic resistance
  !> @param[in]     canhc_tile           Canopy heat capacity on tiles
  !> @param[in]     tile_water_extract   Extraction of water from each tile
  !> @param[in]     ustar                Friction velocity
  !> @param[in,out] lake_evap            Lake evaporation (grid box mean)
  !> @param[in]     soil_moist_avail     Available soil moisture for evaporation
  !> @param[in]     bl_type_ind          Diagnosed BL types
  !> @param[in,out] surf_ht_flux         Surface to sub-surface heat flux
  !> @param[in,out] t1p5m_surft          Diagnostic: 1.5m temperature for land tiles
  !> @param[in,out] q1p5m_surft          Diagnostic: 1.5m specific humidity for land tiles
  !> @param[in,out] t1p5m                Diagnostic: 1.5m temperature
  !> @param[in,out] q1p5m                Diagnostic: 1.5m specific humidity
  !> @param[in,out] qcl1p5m              Diagnostic: 1.5m specific cloud water
  !> @param[in,out] rh1p5m               Diagnostic: 1.5m relative humidity
  !> @param[in,out] t1p5m_ssi            Diagnostic: 1.5m temperature over sea and sea-ice
  !> @param[in,out] q1p5m_ssi            Diagnostic: 1.5m specific humidity over sea and sea-ice
  !> @param[in,out] qcl1p5m_ssi          Diagnostic: 1.5m specific cloud water over sea and sea-ice
  !> @param[in,out] rh1p5m_ssi           Diagnostic: 1.5m relative humidity over sea and sea-ice
  !> @param[in,out] t1p5m_land           Diagnostic: 1.5m temperature over land
  !> @param[in,out] q1p5m_land           Diagnostic: 1.5m specific humidity over land
  !> @param[in,out] qcl1p5m_land         Diagnostic: 1.5m specific cloud water over land
  !> @param[in,out] rh1p5m_land          Diagnostic: 1.5m relative humidity over land
  !> @param[in,out] latent_heat          Diagnostic: Surface latent heat flux
  !> @param[in,out] snomlt_surf_htf      Diagnostic: Grid mean suface snowmelt heat flux
  !> @param[in,out] soil_evap            Diagnostic: Grid mean evapotranspiration from the soil
  !> @param[in,out] soil_surf_ht_flux    Diagnostic: Grid mean surface soil heat flux
  !> @param[in,out] surf_sw_net          Diagnostic: Net surface shortwave radiation
  !> @param[in,out] surf_radnet          Diagnostic: Net surface radiation
  !> @param[in,out] surf_lw_up           Diagnostic: Upward surface longtwave radiation
  !> @param[in,out] surf_lw_down         Diagnostic: Downward surface longwave radiation
  !> @param[in,out] ocn_cpl_point        Diagnostic: Coupling point mask
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_tile             Number of DOFs per cell for tiles
  !> @param[in]     undf_tile            Number of total DOFs for tiles
  !> @param[in]     map_tile             Dofmap for cell for surface tiles
  !> @param[in]     ndf_sice             Number of DOFs per cell for sice levels
  !> @param[in]     undf_sice            Number of total DOFs for sice levels
  !> @param[in]     map_sice             Dofmap for cell for sice levels
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_soil             Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil            Number of total DOFs for soil levels
  !> @param[in]     map_soil             Dofmap for cell for soil levels
  !> @param[in]     ndf_smtile           Number of DOFs per cell for soil levels and tiles
  !> @param[in]     undf_smtile          Number of total DOFs for soil levels and tiles
  !> @param[in]     map_smtile           Dofmap for cell for soil levels and tiles
  !> @param[in]     ndf_bl               Number of DOFs per cell for BL types
  !> @param[in]     undf_bl              Number of total DOFs for BL types
  !> @param[in]     map_bl               Dofmap for cell for BL types
  subroutine jules_imp_code(nlayers, seg_len,                   &
                            outer,                              &
                            loop,                               &
                            wetrho_in_w3,                       &
                            exner_in_wth,                       &
                            height_wth,                         &
                            tile_fraction,                      &
                            sea_ice_thickness,                  &
                            sea_ice_temperature,                &
                            sea_ice_pensolar,                   &
                            sea_ice_pensolar_frac_direct,       &
                            sea_ice_pensolar_frac_diffuse,      &
                            tile_temperature,                   &
                            screen_temperature,                 &
                            time_since_transition,              &
                            latitude,                           &
                            tile_snow_mass,                     &
                            n_snow_layers,                      &
                            snow_depth,                         &
                            canopy_water,                       &
                            soil_temperature,                   &
                            tile_heat_flux,                     &
                            tile_moisture_flux,                 &
                            sw_up_tile,                         &
                            sw_down_surf,                       &
                            lw_down_surf,                       &
                            sw_down_blue_surf,                  &
                            sw_direct_blue_surf,                &
                            skyview,                            &
                            tile_lw_grey_albedo,                &
                            snowice_sublimation,                &
                            surf_heat_flux,                     &
                            canopy_evap,                        &
                            water_extraction,                   &
                            snowice_melt,                       &
                            m_cf,                               &
                            rh_crit_wth,                        &
                            rhokh_bl,                           &
                            moist_flux_bl,                      &
                            heat_flux_bl,                       &
                            dtrdz_tq_bl,                        &
                            alpha1_tile,                        &
                            ashtf_prime_tile,                   &
                            dtstar_tile,                        &
                            fracaero_t_tile,                    &
                            fracaero_s_tile,                    &
                            z0h_tile,                           &
                            z0m_tile,                           &
                            rhokh_tile,                         &
                            chr1p5m_tile,                       &
                            resfs_tile,                         &
                            canhc_tile,                         &
                            tile_water_extract,                 &
                            ustar,                              &
                            lake_evap,                          &
                            soil_moist_avail,                   &
                            bl_type_ind,                        &
                            qw_wth, tl_wth,                     &
                            dqw1_2d, dtl1_2d,                   &
                            ct_ctq1_2d, surf_ht_flux,           &
                            t1p5m_surft, q1p5m_surft,           &
                            t1p5m, q1p5m, qcl1p5m, rh1p5m,      &
                            t1p5m_ssi, q1p5m_ssi,               &
                            qcl1p5m_ssi, rh1p5m_ssi,            &
                            t1p5m_land, q1p5m_land,             &
                            qcl1p5m_land, rh1p5m_land,          &
                            latent_heat, snomlt_surf_htf,       &
                            soil_evap,                          &
                            soil_surf_ht_flux,                  &
                            surf_sw_net, surf_radnet,           &
                            surf_lw_up, surf_lw_down,           &
                            ocn_cpl_point,                      &
                            ndf_w3,                             &
                            undf_w3,                            &
                            map_w3,                             &
                            ndf_wth,                            &
                            undf_wth,                           &
                            map_wth,                            &
                            ndf_tile, undf_tile, map_tile,      &
                            ndf_sice, undf_sice, map_sice,      &
                            ndf_2d,                             &
                            undf_2d,                            &
                            map_2d,                             &
                            ndf_soil, undf_soil, map_soil,      &
                            ndf_smtile, undf_smtile, map_smtile,&
                            ndf_bl, undf_bl, map_bl)

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use gas_calc_all_mod, only: co2_mix_ratio_now
    use jules_control_init_mod, only: n_land_tile, n_sea_ice_tile, &
                                      first_sea_tile, first_sea_ice_tile

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use ancil_info, only: nsurft, nsoilt, dim_cslayer, rad_nband, nmasst
    use atm_fields_bounds_mod, only: pdims, pdims_s
    use atm_step_local, only: dim_cs1
    use bl_option_mod, only: l_noice_in_turb, alpha_cd, flux_bc_opt,          &
                             interactive_fluxes, specified_fluxes_only, puns, &
                             pstb, on
    use carbon_options_mod, only: l_co2_interactive
    use csigma, only: sbcon
    use dust_parameters_mod, only: ndiv, ndivh
    use free_tracers_inputs_mod, only: l_wtrac
    use gen_phys_inputs_mod, only: l_mr_physics
    use jules_deposition_mod, only: l_deposition
    use jules_irrig_mod, only: irr_crop, irr_crop_doell
    use jules_sea_seaice_mod, only: nice, nice_use, emis_sea, l_sice_swpen
    use jules_snow_mod, only: cansnowtile, rho_snow_const, l_snowdep_surf,nsmax
    use jules_surface_types_mod, only: npft, ntype, lake, nnvg, ncpft, nnpft, &
                                       soil
    use jules_surface_mod, only: l_flake_model
    use jules_vegetation_mod, only: can_model, l_crop, l_triffid, l_phenol,    &
                                    can_rad_mod, l_acclim, l_sugar, l_red
    use jules_radiation_mod, only: l_albedo_obs
    use jules_soil_mod, only: ns_deep, l_bedrock
    use jules_soil_biogeochem_mod, only: dim_ch4layer, soil_bgc_model,         &
                                         soil_model_ecosse, l_layeredc
    use jules_water_tracers_mod, only: l_wtrac_jls, n_wtrac_jls, n_evap_srce
    use nlsizes_namelist_mod, only: sm_levels, ntiles, bl_levels
    use planet_constants_mod, only: p_zero, kappa, planet_radius, two_omega
    use rad_input_mod, only: co2_mmr
    use theta_field_sizes, only: t_i_length, t_j_length, &
                                 u_i_length,u_j_length,  &
                                 v_i_length,v_j_length

    ! subroutines used
    use sf_diags_mod, only: dealloc_sf_imp, alloc_sf_imp, strnewsfdiag
    use surf_couple_implicit_mod, only: surf_couple_implicit
    use tilepts_mod, only: tilepts
    use ls_cld_mod, only: ls_cld
    use qsat_mod, only: qsat

    !---------------------------------------
    ! JULES modules
    !---------------------------------------
    use crop_vars_mod,            only: crop_vars_type, crop_vars_data_type,   &
                                        crop_vars_alloc, crop_vars_assoc, &
                                        crop_vars_nullify, crop_vars_dealloc
    use prognostics,              only: progs_data_type, progs_type,           &
                                        prognostics_alloc, prognostics_assoc,  &
                                        prognostics_nullify, prognostics_dealloc
    use jules_vars_mod,           only: jules_vars_type, jules_vars_data_type, &
                                        jules_vars_alloc, jules_vars_assoc,    &
                                        jules_vars_dealloc, jules_vars_nullify
    use aero,                     only: aero_type, aero_data_type,             &
                                        aero_alloc, aero_assoc, &
                                        aero_nullify, aero_dealloc
    use coastal,                  only: coastal_type, coastal_data_type,       &
                                        coastal_assoc, coastal_alloc, &
                                        coastal_dealloc, coastal_nullify
    use lake_mod,                 only: lake_type, lake_data_type,             &
                                        lake_assoc, lake_alloc, &
                                        lake_dealloc, lake_nullify
    use jules_forcing_mod,        only: forcing_type, forcing_data_type,       &
                                        forcing_assoc, forcing_alloc,          &
                                        forcing_nullify, forcing_dealloc
    use ancil_info,               only: ainfo_type, ainfo_data_type,           &
                                        ancil_info_assoc, ancil_info_alloc,    &
                                        ancil_info_dealloc, ancil_info_nullify
    use fluxes_mod,               only: fluxes_type, fluxes_data_type,         &
                                        fluxes_alloc, fluxes_assoc,            &
                                        fluxes_nullify, fluxes_dealloc
    use jules_wtrac_type_mod,     only: jls_wtrac_type, jls_wtrac_data_type,   &
                                        wtrac_jls_assoc, wtrac_jls_alloc,      &
                                        wtrac_jls_nullify, wtrac_jls_dealloc
    use trifctl,                  only: trifctl_type
    ! In general CABLE utilizes a required subset of tbe JULES types, however;
    use progs_cbl_vars_mod, only: progs_cbl_vars_type ! CABLE requires extra progs
    use work_vars_mod_cbl,  only: work_vars_type      ! and some kept thru timestep

    use timestepping_config_mod,       only : method, method_jules

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: outer, loop

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth,seg_len)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3,seg_len)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d,seg_len)

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile,seg_len)
    integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
    integer(kind=i_def), intent(in) :: map_sice(ndf_sice,seg_len)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil,seg_len)
    integer(kind=i_def), intent(in) :: ndf_smtile, undf_smtile
    integer(kind=i_def), intent(in) :: map_smtile(ndf_smtile,seg_len)

    integer(kind=i_def), intent(in) :: ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl,seg_len)

    real(kind=r_def), dimension(undf_w3),  intent(inout) :: moist_flux_bl,     &
                                                            heat_flux_bl
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: wetrho_in_w3,       &
                                                           rhokh_bl
    real(kind=r_def), dimension(undf_wth), intent(in)   :: exner_in_wth,       &
                                                           height_wth,         &
                                                           rh_crit_wth,        &
                                                           dtrdz_tq_bl,        &
                                                           m_cf, qw_wth, tl_wth

    real(kind=r_def), dimension(undf_2d), intent(in) :: ustar,                &
                                                        soil_moist_avail

    real(kind=r_def), intent(in)    :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(inout) :: screen_temperature(undf_tile)
    real(kind=r_def), intent(inout) :: lake_evap(undf_2d)
    real(kind=r_def), intent(inout) :: time_since_transition(undf_2d)
    real(kind=r_def), intent(in)    :: latitude(undf_2d)
    real(kind=r_def), intent(in)    :: tile_snow_mass(undf_tile)
    integer(kind=i_def), intent(in) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(in)    :: snow_depth(undf_tile)
    real(kind=r_def), intent(in)    :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_heat_flux(undf_tile)
    real(kind=r_def), intent(inout) :: tile_moisture_flux(undf_tile)
    real(kind=r_def), intent(in)    :: sw_up_tile(undf_tile)
    real(kind=r_def), intent(inout)   :: snowice_sublimation(undf_tile)
    real(kind=r_def), intent(inout)   :: surf_heat_flux(undf_tile)
    real(kind=r_def), intent(inout)   :: canopy_evap(undf_tile)
    real(kind=r_def), intent(inout)   :: snowice_melt(undf_tile)

    real(kind=r_def), intent(in)    :: sea_ice_thickness(undf_sice)
    real(kind=r_def), intent(inout) :: sea_ice_temperature(undf_sice)
    real(kind=r_def), intent(inout) :: sea_ice_pensolar(undf_sice)
    real(kind=r_def), intent(in)    :: sea_ice_pensolar_frac_direct(undf_sice)
    real(kind=r_def), intent(in)    :: sea_ice_pensolar_frac_diffuse(undf_sice)

    real(kind=r_def), intent(in) :: sw_down_surf(undf_2d)
    integer(kind=i_def), intent(in) :: ocn_cpl_point(undf_2d)
    real(kind=r_def), intent(in) :: lw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_blue_surf(undf_2d)
    real(kind=r_def), intent(in) :: sw_direct_blue_surf(undf_2d)
    real(kind=r_def), intent(in) :: skyview(undf_2d)
    real(kind=r_def), intent(in) :: tile_lw_grey_albedo(undf_tile)
    real(kind=r_def), intent(in) :: dqw1_2d(undf_2d)
    real(kind=r_def), intent(in) :: dtl1_2d(undf_2d)
    real(kind=r_def), intent(in) :: ct_ctq1_2d(undf_2d)
    real(kind=r_def), intent(inout) :: surf_ht_flux(undf_tile)
    real(kind=r_def), pointer, intent(inout) :: t1p5m_surft(:)
    real(kind=r_def), pointer, intent(inout) :: q1p5m_surft(:)
    real(kind=r_def), pointer, intent(inout) :: t1p5m(:), q1p5m(:)
    real(kind=r_def), pointer, intent(inout) :: qcl1p5m(:), rh1p5m(:)
    real(kind=r_def), pointer, intent(inout) :: t1p5m_ssi(:), q1p5m_ssi(:)
    real(kind=r_def), pointer, intent(inout) :: qcl1p5m_ssi(:), rh1p5m_ssi(:)
    real(kind=r_def), pointer, intent(inout) :: t1p5m_land(:), q1p5m_land(:)
    real(kind=r_def), pointer, intent(inout) :: qcl1p5m_land(:), rh1p5m_land(:)
    real(kind=r_def), pointer, intent(inout) :: latent_heat(:)
    real(kind=r_def), pointer, intent(inout) :: snomlt_surf_htf(:)
    real(kind=r_def), pointer, intent(inout) :: soil_evap(:)
    real(kind=r_def), pointer, intent(inout) :: soil_surf_ht_flux(:)
    real(kind=r_def), pointer, intent(inout) :: surf_sw_net(:)
    real(kind=r_def), pointer, intent(inout) :: surf_radnet(:)
    real(kind=r_def), pointer, intent(inout) :: surf_lw_up(:)
    real(kind=r_def), pointer, intent(inout) :: surf_lw_down(:)

    real(kind=r_def), intent(in) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(inout):: water_extraction(undf_soil)

    real(kind=r_def), intent(in) :: tile_water_extract(undf_smtile)

    integer(kind=i_def), dimension(undf_bl), intent(in) :: bl_type_ind
    real(kind=r_def), dimension(undf_tile), intent(in)  :: alpha1_tile,      &
                                                           ashtf_prime_tile, &
                                                           dtstar_tile,      &
                                                           fracaero_t_tile,  &
                                                           fracaero_s_tile,  &
                                                           z0h_tile,         &
                                                           z0m_tile,         &
                                                           rhokh_tile,       &
                                                           chr1p5m_tile,     &
                                                           resfs_tile,       &
                                                           canhc_tile

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: i, i_tile, i_sice, n, l, m, land_field, ssi_pts,       &
                      sice_pts, sea_pts, k

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    integer(i_um) :: river_row_length_dum, river_rows_dum

    ! local switches and scalars
    integer(i_um) :: error_code
    logical :: l_correct
    logical :: l_wtrac_bl

    real(r_def) :: sw_diffuse_blue_surf

    ! profile fields from level 1 upwards
    real(r_um), dimension(seg_len,1) :: rhcpt, qcf_latest, co2

    ! profile field on boundary layer levels
    real(r_um), dimension(seg_len,1) :: fqw, ftl, rhokh, rhokh_mix

    ! profile fields on u/v points and all levels
    real(r_um), dimension(seg_len,1) :: u, v

    ! profile fields on u/v points and BL levels
    real(r_um), dimension(seg_len,1) :: taux, tauy, rhokm_u, rhokm_v, &
         du_star, dv_star

    ! single level real fields
    real(r_um), dimension(seg_len,1) :: tstar_land, dtstar_sea,              &
         tstar_sice, alpha1_sea, ashtf_prime_sea, chr1p5m_sice, flandg,      &
         rhokh_sea, u_s, z0hssi, z0mssi, work_2d_1, work_2d_2, work_2d_3,    &
         qcl1p5m_loc, gamma1, gamma2, ctctq1_1, dqw1_1, dtl1_1, cq_cm_u_1,   &
         du_1, cq_cm_v_1, dv_1, tscrndcl_ssi, tstbtrans, f3_at_p, rho1,      &
         ti_sice, olr, sky

    ! single level real fields on u/v points
    real(r_um), dimension(seg_len,1) :: flandg_u, flandg_v, cdr10m_u, cdr10m_v

    ! single level integer fields
    integer(i_um), dimension(seg_len,1) :: ntml

    ! single level logical fields
    logical, dimension(seg_len,1) :: cumulus

    ! fields on sea-ice categories
    real(r_um), dimension(seg_len,1,nice_use) :: alpha1_sice, &
         ashtf_prime, rhokh_sice, dtstar_sice, radnet_sice

    ! fields on land points and surface tiles
    real(r_um), dimension(:,:), allocatable :: dtstar_surft,              &
         alpha1, ashtf_prime_surft, chr1p5m,                              &
         fracaero_t, fracaero_s, resfs, rhokh_surft,                      &
         resft, flake, canhc_surft, tscrndcl_surft, epot_surft

    ! field on surface tiles and soil levels
    real(r_um), dimension(:,:,:), allocatable :: wt_ext_surft

    ! fields on all points
    real(r_um), dimension(:,:), allocatable :: t1p5m_land_loc, q1p5m_land_loc

    ! parameters for new BL solver
    real(r_um) :: pnonl,p1,p2
    real(r_um), dimension(seg_len) :: i1, e1, e2
    real(r_um), parameter :: sqrt2 = sqrt(2.0_r_def)
    real(r_um) :: r_gamma

    !-----------------------------------------------------------------------
    ! JULES Types
    !-----------------------------------------------------------------------
    type(strnewsfdiag)         :: sf_diag
    type(crop_vars_type)       :: crop_vars
    type(crop_vars_data_type)  :: crop_vars_data
    type(progs_type)           :: progs
    type(progs_data_type)      :: progs_data
    type(jules_vars_type)      :: jules_vars
    type(jules_vars_data_type) :: jules_vars_data
    type(aero_type)            :: aerotype
    type(aero_data_type)       :: aero_data
    type(coastal_type)         :: coast
    type(coastal_data_type)    :: coastal_data
    type(lake_type)            :: lake_vars
    type(lake_data_type)       :: lake_data
    type(ainfo_type)           :: ainfo
    type(ainfo_data_type)      :: ainfo_data
    type(forcing_type)         :: forcing
    type(forcing_data_type)    :: forcing_data
    type(fluxes_type)          :: fluxes
    type(fluxes_data_type)     :: fluxes_data
    type(jls_wtrac_type)       :: wtrac_jls
    type(jls_wtrac_data_type)  :: wtrac_jls_data
    !CABLE TYPES containing field data (IN OUT)
    type(progs_cbl_vars_type) :: progs_cbl_vars
    type(work_vars_type)      :: work_cbl


    !-----------------------------------------------------------------------
    ! Initialisation of JULES data and pointer types
    !-----------------------------------------------------------------------

    ! Land tile fractions
    land_field = 0
    do i = 1, seg_len
      flandg(i,1) = 0.0_r_um
      do n = 1, n_land_tile
        flandg(i,1) = flandg(i,1) + real(tile_fraction(map_tile(1,i)+n-1), r_um)
      end do
      if (flandg(i,1) > 0.0_r_um) then
        land_field = land_field + 1
      end if
    end do

    call crop_vars_alloc(land_field, t_i_length, t_j_length,                  &
                     nsurft, ncpft,nsoilt, sm_levels, l_crop, irr_crop,       &
                     irr_crop_doell, crop_vars_data)

    call crop_vars_assoc(crop_vars, crop_vars_data)

    call prognostics_alloc(land_field, t_i_length, t_j_length,                &
                      nsurft, npft, nsoilt, sm_levels, ns_deep, nsmax,        &
                      dim_cslayer, dim_cs1, dim_ch4layer,                     &
                      nice, nice_use, soil_bgc_model, soil_model_ecosse,      &
                      l_layeredc, l_triffid, l_phenol, l_bedrock, l_red,      &
                      nmasst, nnpft, l_acclim, l_sugar, progs_data)
    call prognostics_assoc(progs,progs_data)

    call jules_vars_alloc(land_field,ntype,nsurft,rad_nband,nsoilt,sm_levels, &
                t_i_length, t_j_length, npft, bl_levels, pdims_s, pdims,      &
                l_albedo_obs, cansnowtile, l_deposition,                      &
                jules_vars_data)
    call jules_vars_assoc(jules_vars,jules_vars_data)

    call aero_alloc(land_field,t_i_length,t_j_length,                         &
                nsurft,ndiv, aero_data)
    call aero_assoc(aerotype, aero_data)

    call coastal_alloc(land_field,t_i_length,t_j_length,                      &
                   u_i_length,u_j_length,                                     &
                   v_i_length,v_j_length,                                     &
                   nice_use,nice,coastal_data)
    call coastal_assoc(coast, coastal_data)

    call lake_alloc(land_field, l_flake_model, lake_data)
    call lake_assoc(lake_vars, lake_data)

    call ancil_info_alloc(land_field,t_i_length,t_j_length,                   &
                      nice,nsoilt,ntype,                                      &
                      ainfo_data)
    call ancil_info_assoc(ainfo, ainfo_data)

    call forcing_alloc(t_i_length,t_j_length,u_i_length, u_j_length,          &
                       v_i_length, v_j_length, forcing_data)
    call forcing_assoc(forcing, forcing_data)

    call fluxes_alloc(land_field, t_i_length, t_j_length,                     &
                      nsurft, npft, nsoilt, sm_levels,                        &
                      nice, nice_use,                                         &
                      fluxes_data)
    call fluxes_assoc(fluxes, fluxes_data)

    ! Set river size to 1 here as the fields with these dimensions are not
    ! used here
    river_row_length_dum = 1
    river_rows_dum = 1
    call wtrac_jls_alloc(land_field, seg_len, 1, n_land_tile, nsoilt,          &
                         sm_levels, nsmax, nice_use, n_wtrac_jls, n_evap_srce, &
                         river_row_length_dum, river_rows_dum, l_wtrac_jls,    &
                         wtrac_jls_data)
    call wtrac_jls_assoc(wtrac_jls, wtrac_jls_data)

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    !-----------------------------------------------------------------------
    error_code=0

    ! Set logical flags for sf_diags
    sf_diag%smlt = .not. associated(snomlt_surf_htf, empty_real_data)
    sf_diag%l_lw_surft = .not. associated(surf_lw_up, empty_real_data)  &
                    .or. .not. associated(surf_lw_down, empty_real_data)
    sf_diag%l_lw_up_sice_weighted_cat = .not. associated(surf_lw_up, empty_real_data)
    sf_diag%sq1p5 = .true.
    call alloc_sf_imp(sf_diag, outer == outer_iterations, land_field)
    if (sf_diag%simlt) then
      do n = 1, nice
        do i = 1, seg_len
          sf_diag%sice_mlt_htf(i,1,n) = 0.0_r_um
        end do
      end do
    end if

    allocate(epot_surft(land_field,ntiles))
    allocate(tscrndcl_surft(land_field,ntiles))
    allocate(alpha1(land_field,ntiles))
    allocate(ashtf_prime_surft(land_field,ntiles))
    allocate(dtstar_surft(land_field,ntiles))
    allocate(fracaero_t(land_field,ntiles))
    allocate(fracaero_s(land_field,ntiles))
    allocate(rhokh_surft(land_field,ntiles))
    allocate(chr1p5m(land_field,ntiles))
    allocate(resfs(land_field,ntiles))
    allocate(canhc_surft(land_field,ntiles))
    allocate(resft(land_field,ntiles))
    allocate(flake(land_field,ntiles))
    allocate(wt_ext_surft(land_field,sm_levels,ntiles))

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------

    ! Sea-ice fraction
    do i = 1, seg_len
      do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = n - first_sea_ice_tile + 1
        ainfo%ice_fract_ij(i,1) = ainfo%ice_fract_ij(i,1) + &
             real(tile_fraction(map_tile(1,i)+n-1), r_um)
        ainfo%ice_fract_ncat_sicat(i, 1, i_sice) = real(tile_fraction(map_tile(1,i)+n-1), r_um)
      end do

      ! Because Jules tests on flandg < 1, we need to ensure this is exactly
      ! 1 when no sea or sea-ice is present
      if ( tile_fraction(map_tile(1,i)+first_sea_tile-1) == 0.0_r_def .and. &
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

    do i = 1, seg_len
      fqw(i,1) = moist_flux_bl(map_w3(1,i))
      ftl(i,1) = heat_flux_bl(map_w3(1,i))
    end do

    if (loop == 1) then

      l = 0
      do i = 1, seg_len
        if (flandg(i,1) > 0.0_r_um) then
          l = l+1
          ainfo%land_index(l) = i
        end if
      end do

      do l = 1, land_field
        coast%fland(l) = flandg(ainfo%land_index(l),1)
        do n = 1, n_land_tile
          ! Jules requires fractions with respect to the land area
          ainfo%frac_surft(l, n) = real(tile_fraction(map_tile(1,ainfo%land_index(l))+n-1), r_um) &
             / coast%fland(l)
        end do
      end do

      ! Set type_pts and type_index
      call tilepts(land_field, ainfo%frac_surft, ainfo%surft_pts,              &
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
      sice_pts = 0
      sea_pts = 0
      ! then calculate based on state
      do i = 1, seg_len
        if (ainfo%ssi_index(i) > 0) then
          if (ainfo%ice_fract_ij(i, 1) > 0.0_r_um) then
            sice_pts = sice_pts + 1
            ainfo%sice_index(sice_pts) = i
            ainfo%sice_frac(i) = ainfo%ice_fract_ij(i, 1)
          end if
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
          tscrndcl_surft(l, n) = real(screen_temperature(map_tile(1,ainfo%land_index(l))+n-1), r_um)
          tstar_land(ainfo%land_index(l),1) = tstar_land(ainfo%land_index(l),1) &
             + ainfo%frac_surft(l, n) * progs%tstar_surft(l, n)
          ! sensible heat flux
          fluxes%ftl_surft(l, n) = real(tile_heat_flux(map_tile(1,ainfo%land_index(l))+n-1), r_um)
          ! moisture flux
          fluxes%fqw_surft(l, n) = real(tile_moisture_flux(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        end do
      end do

      ! Fields for decoupled screen temperature diagnostic
      do i = 1, seg_len
        tstbtrans(i,1)  = time_since_transition(map_2d(1,i))
        f3_at_p(i,1) = two_omega * sin(latitude(map_2d(1,i)))
      end do

      ! Sea temperature
      ! Default to temperature over frozen sea as the initialisation
      ! that follows does not initialise sea points if they are fully
      ! frozen
      do i = 1, seg_len
        if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then
          coast%tstar_sea_ij(i,1) = real(tile_temperature(map_tile(1,i)+first_sea_tile-1), r_um)
        else
          coast%tstar_sea_ij(i,1) = tfs
        end if
        tscrndcl_ssi(i,1)=real(screen_temperature(map_tile(1,i)+first_sea_tile-1), r_um)
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
            ! sea-ice heat flux
            fluxes%ftl_sicat(i,1,i_sice) = real(tile_heat_flux(map_tile(1,i)+n-1), r_um)
            ! sea-ice moisture flux
            fluxes%fqw_sicat(i,1,i_sice) = real(tile_moisture_flux(map_tile(1,i)+n-1), r_um)
          end do
        end if
      end do

      do i = 1, seg_len
        ! Sea & Sea-ice temperature
        coast%tstar_ssi_ij(i,1) = (1.0_r_um - ainfo%ice_fract_ij(i,1)) &
           * coast%tstar_sea_ij(i,1) + ainfo%ice_fract_ij(i,1) * tstar_sice(i,1)
      end do

      if (can_rad_mod == 6) then
        jules_vars%diff_frac = 0.4_r_um
      end if

      ! Following variables need to be initialised to stop crashed in unused
      ! UM code
      do i = 1, seg_len
        olr(i,1) = 300.0_r_um
        cdr10m_u(i,1) = 0.0_r_um
        cdr10m_v(i,1) = 0.0_r_um
      end do
      do n = 1, ntiles
        do l = 1, land_field
          epot_surft(l,n) = 0.0_r_um
        end do
      end do

      do l = 1, land_field
        do m = 1, sm_levels
          ! Soil temperature (t_soil_soilt)
          progs%t_soil_soilt(l, 1, m) = real(soil_temperature(map_soil(1,ainfo%land_index(l))+m-1), r_um)
        end do
      end do

      if (topography == topography_horizon) then
        ! Set skyview factor used internally by JULES
        do i = 1, seg_len
          sky(i,1) = real(skyview(map_2d(1,i)), r_um)
        end do
      end if

      do i = 1, seg_len
        ! Downwelling LW radiation at surface
        forcing%lw_down_ij(i,1) = real(lw_down_surf(map_2d(1,i)), r_um)
      end do

      ! Net SW radiation on tiles
      do l = 1, land_field
        do n = 1, n_land_tile
          fluxes%sw_surft(l, n) = real(sw_down_surf(map_2d(1,ainfo%land_index(l))) - &
                                sw_up_tile(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        end do
      end do

      do i = 1, seg_len

        ! The amount of diffuse visible light is needed for solar penetrating radiation
        sw_diffuse_blue_surf = sw_down_blue_surf(map_2d(1,i)) - sw_direct_blue_surf(map_2d(1,i))

        ! Calculate penetrating solar into sea ice. We use the fraction of penetrating
        ! solar for direct visible light as it is the same as for diffuse visible light.
        if (l_sice_swpen .and. ocn_cpl_point(map_2d(1,i)) == 1_i_def ) then
          do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
            i_sice = n - first_sea_ice_tile + 1
            sea_ice_pensolar(map_sice(1,i)+i_sice-1) = sw_direct_blue_surf(map_2d(1,i)) *    &
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
      end do

      ! Carbon dioxide
      co2_mmr = real(co2_mix_ratio_now, r_um)
      do i = 1, seg_len
        co2(i,1) = co2_mmr
      end do

      do l = 1, land_field
        ! Canopy water on each tile (canopy_surft)
        do n = 1, n_land_tile
          progs%canopy_surft(l, n) = real(canopy_water(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        end do

        do n = 1, n_land_tile
          ! Lying snow mass on land tiles
          progs%snow_surft(l, n) = real(tile_snow_mass(map_tile(1,ainfo%land_index(l))+n-1), r_um)
          ! Number of snow layers on tiles (nsnow_surft)
          progs%nsnow_surft(l, n) = n_snow_layers(map_tile(1,ainfo%land_index(l))+n-1)
        end do
      end do
      ! Equivalent snowdepth for surface calculations.
      ! code copied from jules_land_sf_explicit
      ! 4 is a magic number inherited from Jules, meaning radiative canopy
      ! with heat capacity and snow beneath
      do n = 1, n_land_tile
        if ( (can_model == 4) .and. cansnowtile(n) .and. l_snowdep_surf) then
          do l = 1, land_field
            jules_vars%snowdep_surft(l, n) = progs%snow_surft(l, n) / rho_snow_const
          end do
        else
          do l = 1, land_field
            jules_vars%snowdep_surft(l, n) = real(snow_depth(map_tile(1,ainfo%land_index(l))+n-1), r_um)
          end do
        end if
      end do

      !-----------------------------------------------------------------------
      ! Assuming map_wth(1) points to level 0
      ! and map_w3(1) points to level 1
      !-----------------------------------------------------------------------
      do i = 1, seg_len
        ! Density on lowest level
        rho1(i,1) = wetrho_in_w3(map_w3(1,i))
        ! surface pressure
        forcing%pstar_ij(i,1) = p_zero*(exner_in_wth(map_wth(1,i) + 0))**(1.0_r_def/kappa)
        ! height of lowest level above surface
        ainfo%z1_tq_ij(i,1) = height_wth(map_wth(1,i) + 1) &
                            - height_wth(map_wth(1,i) + 0)
      end do

      !-----------------------------------------------------------------------
      ! Things passed from other parametrization schemes on this timestep
      !-----------------------------------------------------------------------
      do i = 1, seg_len
        rhokh(i,1) = rhokh_bl(map_w3(1,i))
        forcing%qw_1_ij(i,1) = qw_wth(map_wth(1,i) + 1)
        forcing%tl_1_ij(i,1) = tl_wth(map_wth(1,i) + 1)
        jules_vars%dtrdz_charney_grid_1_ij(i,1) = dtrdz_tq_bl(map_wth(1,i) + 1)
      end do

      do l = 1, land_field
        do n = 1, n_land_tile
          alpha1(l, n) = alpha1_tile(map_tile(1,ainfo%land_index(l))+n-1)
          fracaero_t(l, n) = fracaero_t_tile(map_tile(1,ainfo%land_index(l))+n-1)
          fracaero_s(l, n) = fracaero_s_tile(map_tile(1,ainfo%land_index(l))+n-1)
          fluxes%z0h_surft(l, n) = z0h_tile(map_tile(1,ainfo%land_index(l))+n-1)
          fluxes%z0m_surft(l, n) = z0m_tile(map_tile(1,ainfo%land_index(l))+n-1)
          chr1p5m(l, n) = chr1p5m_tile(map_tile(1,ainfo%land_index(l))+n-1)
          resfs(l, n) = resfs_tile(map_tile(1,ainfo%land_index(l))+n-1)
          canhc_surft(l, n) = canhc_tile(map_tile(1,ainfo%land_index(l))+n-1)
          ! recalculate the surface emissivity
          fluxes%emis_surft(l, n) = 1.0_r_um - &
               real(tile_lw_grey_albedo(map_tile(1,ainfo%land_index(l))+n-1), r_um)
        end do
        ! recalculate the total resistance factor
        do n = 1, n_land_tile
          resft(l,n) = fracaero_t(l,n) + (1.0_r_um - fracaero_t(l,n)) * resfs(l,n)
        end do
        resft(l,lake) = 1.0_r_um
        ! recalculate the total lake fraction
        do n = 1, n_land_tile
          flake(l,n) = 0.0_r_um
        end do
        flake(l,lake) = 1.0_r_um
      end do
      ! Fields only calculated on tiles which exist
      do n = 1, n_land_tile
        do m = 1, ainfo%surft_pts(n)
          l = ainfo%surft_index(m,n)
          ashtf_prime_surft(l, n) = ashtf_prime_tile(map_tile(1,ainfo%land_index(l))+n-1)
          dtstar_surft(l, n) = dtstar_tile(map_tile(1,ainfo%land_index(l))+n-1)
          rhokh_surft(l, n) = rhokh_tile(map_tile(1,ainfo%land_index(l))+n-1)
        end do
      end do

      do n = 1, n_land_tile
        do m = 1, sm_levels
          i_tile=sm_levels*(n-1) + m - 1  
          do k = 1, ainfo%surft_pts(n)
            l = ainfo%surft_index(k,n)
            wt_ext_surft(l,m,n) = tile_water_extract(map_smtile(1,ainfo%land_index(l))+i_tile)
          end do
        end do
      end do

      do i = 1, seg_len
        alpha1_sea(i,1) = alpha1_tile(map_tile(1,i)+first_sea_tile-1)
      end do
      ! This is only calculated on sea points
      do l = 1, sea_pts
        i = ainfo%sea_index(l)
        ashtf_prime_sea(i,1) = ashtf_prime_tile(map_tile(1,i)+first_sea_tile-1)
        dtstar_sea(i,1) = dtstar_tile(map_tile(1,i)+first_sea_tile-1)
      end do
      ! This is only used with multi-category sea-ice
      if (n_sea_ice_tile > 1) then
        do i = 1, seg_len
          rhokh_sea(i,1) = rhokh_tile(map_tile(1,i)+first_sea_tile-1)
        end do
      end if

      do i = 1, seg_len
        do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = n - first_sea_ice_tile + 1
          alpha1_sice(i,1,i_sice) = alpha1_tile(map_tile(1,i)+n-1)
          rhokh_sice(i,1,i_sice) = rhokh_tile(map_tile(1,i)+n-1)
        end do
        z0hssi(i,1) = z0h_tile(map_tile(1,i)+first_sea_ice_tile-1)
        z0mssi(i,1) = z0m_tile(map_tile(1,i)+first_sea_ice_tile-1)
        chr1p5m_sice(i,1) = chr1p5m_tile(map_tile(1,i)+first_sea_ice_tile-1)
      end do
      ! Fields are only written on sea-ice points
      do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = n - first_sea_ice_tile + 1
        do l = 1, ainfo%sice_pts_ncat(i_sice)
          i = ainfo%sice_index_ncat(l,i_sice)
          ashtf_prime(i,1,i_sice) = ashtf_prime_tile(map_tile(1,i)+n-1)
          dtstar_sice(i,1,i_sice) = dtstar_tile(map_tile(1,i)+n-1)
        end do
      end do

      do i = 1, seg_len
        u_s(i,1) = ustar(map_2d(1,i))
        dqw1_1(i,1) = dqw1_2d(map_2d(1,i))
        dtl1_1(i,1) = dtl1_2d(map_2d(1,i))
        ctctq1_1(i,1) = ct_ctq1_2d(map_2d(1,i))
      end do

      do l = 1, land_field
        progs%smc_soilt(l,1) = soil_moist_avail(map_2d(1,ainfo%land_index(l)))
      end do

      !-----------------------------------------------------------------------
      ! Fields for diagnostics only
      !-----------------------------------------------------------------------
      if (outer==outer_iterations) then

        if (l_noice_in_turb) then
          do i = 1, seg_len
            qcf_latest(i,1) = 0.0_r_def
          end do
        else
          do i = 1, seg_len
            qcf_latest(i,1) = m_cf(map_wth(1,i) + 1)
          end do
        end if
      end if

    else ! loop=2

      ! Ocean coupling point
      do i = 1, seg_len
        if (ocn_cpl_point(map_2d(1,i)) == 1_i_def) then
          ainfo%ocn_cpl_point(i,1) = .true.
        else
          ainfo%ocn_cpl_point(i,1) = .false.
        end if
      end do

      ! Sea-ice bulk temperature and thickness
      do i = 1, seg_len
        if (ainfo%ice_fract_ij(i, 1) > 0.0_r_um) then
          do n = 1, n_sea_ice_tile
            progs%ti_sicat(i, 1, n) = real(sea_ice_temperature(map_sice(1,i)+n-1), r_um)
            progs%di_ncat_sicat(i, 1, n) = real(sea_ice_thickness(map_sice(1,i)+n-1), r_um)
          end do
        end if
      end do

      if (outer == outer_iterations) then
        do i = 1, seg_len
          do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
            i_sice = n - first_sea_ice_tile + 1
            fluxes%surf_ht_flux_sice(i,1,i_sice) = surf_ht_flux(map_tile(1,i)+n-1)
          end do
        end do
      end if

    end if ! loop

    !-----------------------------------------------------------------------
    ! Code from ni_imp_ctl
    !-----------------------------------------------------------------------

    r_gamma=alpha_cd(1)
    if (flux_bc_opt > interactive_fluxes) then
      ! explicit scalar surface fluxes means surface evaporation, snow melt,
      ! etc will all be consistent with specified surface fluxes
      r_gamma=0.0_r_def
    end if ! flux_bc_opt > interactive_fluxes

    if (method == method_jules) then
      ! Set weights so that coupling is fully explicit,
      ! i.e. gamma2 will be 0 on the 1st call and 1 on the 2nd call
      do i = 1, seg_len
        i1(i) = 1.0_r_def
        e1(i) = 1.0_r_def
        e2(i) = 0.0_r_def
        gamma1(i,1) = r_gamma
      end do
    else
      do i = 1, seg_len
        p1=bl_type_ind(map_bl(1,i)+0)*pstb+(1.0_r_def-bl_type_ind(map_bl(1,i)+0))*puns
        p2=bl_type_ind(map_bl(1,i)+1)*pstb+(1.0_r_def-bl_type_ind(map_bl(1,i)+1))*puns
        pnonl=max(p1,p2)
        i1(i) = (1.0_r_def+1.0_r_def/sqrt2)*(1.0_r_def+pnonl)
        e1(i) = (1.0_r_def+1.0_r_def/sqrt2)*( pnonl + (1.0_r_def/sqrt2) + &
             sqrt(pnonl*(sqrt2-1.0_r_def)+0.5_r_def) )
        e2(i) = (1.0_r_def+1.0_r_def/sqrt2)*( pnonl+(1.0_r_def/sqrt2) - &
             sqrt(pnonl*(sqrt2-1.0_r_def)+0.5_r_def))
        gamma1(i,1) = i1(i)
      end do
    end if

    if (loop == 1) then
      do i = 1, seg_len
        gamma2(i,1) = i1(i) - e1(i)
      end do
      l_correct = .false.
    else
      do i = 1, seg_len
        gamma2(i,1) = i1(i) - e2(i)
      end do
      l_correct = .true.
    end if

    ! Water tracers are only updated on final loop
    l_wtrac_bl = (l_wtrac .AND. outer == outer_iterations)

    call surf_couple_implicit(                                               &
         !Important switch
         l_correct, u, v,                                                    &
         !Misc INTENT(IN) Many of these come out of explicit and into here.
         !Things with _u/v get interpolated in UM
         rhokm_u, rhokm_v, gamma1, gamma2, alpha1, alpha1_sea, alpha1_sice,  &
         ashtf_prime, ashtf_prime_sea, ashtf_prime_surft, du_1, dv_1,        &
         fracaero_t, fracaero_s, resfs, resft,                               &
         rhokh, rhokh_surft, rhokh_sice, rhokh_sea, z0hssi,                  &
         z0mssi, chr1p5m, chr1p5m_sice, canhc_surft, flake, ainfo%frac_surft,&
         wt_ext_surft, cdr10m_u, cdr10m_v, r_gamma,                          &
         !INOUT diagnostics
         sf_diag,                                                            &
         !Fluxes INTENT(INOUT)
         fqw, ftl,                                                           &
         !Misc INTENT(INOUT)
         epot_surft, dtstar_surft, dtstar_sea, dtstar_sice, radnet_sice, olr,&
         !Fluxes INTENT(OUT)
         taux, tauy,                                                         &
         !Misc INTENT(OUT)
         error_code,                                                         &
         !UM-only arguments
         !JULES ancil_info module
         !IN
         ntiles, land_field, ssi_pts, sice_pts, sea_pts, ainfo%surft_pts,    &
         !JULES coastal module IN
         flandg,                                                             &
         ! Coastal OUT - do this here as needs to be OUT
         tstar_land, tstar_sice,                                             &
         ! Water tracer switch (IN)
         l_wtrac_bl,                                                         &
         !JULES switches module
         !IN
         l_co2_interactive, l_mr_physics, co2,                               &
         !Arguments without a JULES module
         !IN
         ctctq1_1,dqw1_1,dtl1_1,du_star,dv_star,cq_cm_u_1,cq_cm_v_1,         &
         flandg_u,flandg_v,                                                  &
         !3 IN, 3 INOUT
         rho1, f3_at_p, u_s,TScrnDcl_SSI,TScrnDcl_SURFT,tStbTrans,           &
         !OUT
         rhokh_mix, ti_sice, sky,                                            &
         !TYPES containing field data (IN OUT)
         crop_vars,ainfo, aerotype, progs, coast, jules_vars,                &
         fluxes, lake_vars, forcing, wtrac_jls, progs_cbl_vars, work_cbl     &
         )

    do i = 1, seg_len
      heat_flux_bl(map_w3(1,i)) = ftl(i,1)
    end do

    if (loop == 1) then
      do i = 1, seg_len
        moist_flux_bl(map_w3(1,i)) = fqw(i,1)
      end do
    end if

    ! Diagnostics
    if (outer == outer_iterations) then

      if (loop==1) then

        do i = 1, seg_len
          tile_heat_flux(map_tile(1,i)+first_sea_tile-1) = 0.0_r_def
          tile_moisture_flux(map_tile(1,i)+first_sea_tile-1) = 0.0_r_def
        end do
        ! Update land tiles
        do n = 1, n_land_tile
          do l = 1, land_field
            tile_temperature(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(progs%tstar_surft(l, n), r_def)
            screen_temperature(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(tscrndcl_surft(l, n), r_def)
            tile_heat_flux(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(fluxes%ftl_surft(l, n), r_def)
            tile_moisture_flux(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(fluxes%fqw_surft(l, n), r_def)
            ! Sum the fluxes over the land for use in sea point calculation
            if (tile_fraction(map_tile(1,ainfo%land_index(l))+n-1) > 0.0_r_def) then
              tile_heat_flux(map_tile(1,ainfo%land_index(l))+first_sea_tile-1) &
                 = tile_heat_flux(map_tile(1,ainfo%land_index(l))+first_sea_tile-1) &
                 + fluxes%ftl_surft(l,n) * tile_fraction(map_tile(1,ainfo%land_index(l))+n-1)
              tile_moisture_flux(map_tile(1,ainfo%land_index(l))+first_sea_tile-1) &
                 = tile_moisture_flux(map_tile(1,ainfo%land_index(l))+first_sea_tile-1) &
                 + fluxes%fqw_surft(l,n) * tile_fraction(map_tile(1,ainfo%land_index(l))+n-1)
            end if
            snowice_sublimation(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(fluxes%ei_surft(l, n), r_def)
            ! NB - net surface heat flux
            surf_heat_flux(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(fluxes%surf_htf_surft(l, n), r_def)
            canopy_evap(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(fluxes%ecan_surft(l, n), r_def)
            snowice_melt(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(fluxes%melt_surft(l, n), r_def)
          end do
        end do

        ! Update sea tile
        do i = 1, seg_len
          tile_temperature(map_tile(1,i)+first_sea_tile-1) = real(coast%tstar_sea_ij(i,1), r_def)
          screen_temperature(map_tile(1,i)+first_sea_tile-1) = real(tscrndcl_ssi(i, 1), r_def)
          time_since_transition(map_2d(1,i)) = tstbtrans(i,1)
        end do

        ! Update sea-ice tiles
        do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = n - first_sea_ice_tile + 1
          do i = 1, seg_len
            tile_temperature(map_tile(1,i)+n-1) = real(coast%tstar_sice_sicat(i,1,i_sice), r_def)
            tile_heat_flux(map_tile(1,i)+n-1) = real(fluxes%ftl_sicat(i,1,i_sice), r_def)
            tile_moisture_flux(map_tile(1,i)+n-1) = real(fluxes%fqw_sicat(i,1,i_sice), r_def)
            snowice_melt(map_tile(1,i)+n-1) = real(fluxes%sice_melt(i,1, i_sice), r_def)
            snowice_sublimation(map_tile(1,i)+n-1) = real(fluxes%ei_sice(i,1,i_sice), r_def)
            ! Sum the fluxes over the sea-ice for use in sea point calculation
            tile_heat_flux(map_tile(1,i)+first_sea_tile-1) =                   &
                 tile_heat_flux(map_tile(1,i)+first_sea_tile-1)                &
                 + tile_heat_flux(map_tile(1,i)+n-1) * tile_fraction(map_tile(1,i)+n-1)
            tile_moisture_flux(map_tile(1,i)+first_sea_tile-1) =               &
                 tile_moisture_flux(map_tile(1,i)+first_sea_tile-1)            &
                 + tile_moisture_flux(map_tile(1,i)+n-1) * tile_fraction(map_tile(1,i)+n-1)
          end do
        end do

        do n = 1, n_land_tile
          do l = 1, land_field
            surf_ht_flux(map_tile(1,ainfo%land_index(l))+n-1) = &
                 real(fluxes%surf_htf_surft(l, n), r_def)
          end do
        end do
        do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = n - first_sea_ice_tile + 1
          do i = 1, seg_len
            if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_def) then
              surf_ht_flux(map_tile(1,i)+n-1) =                             &
                   real(fluxes%surf_ht_flux_sice(i,1,i_sice), r_def)
            end if
          end do
        end do

        do m = 1, sm_levels
          do l = 1, land_field
            water_extraction(map_soil(1,ainfo%land_index(l))+m-1) = &
                 real(fluxes%ext_soilt(l, 1, m), r_def)
          end do
        end do

        ! Lake evaporation is needed for the lake water correction but it is converted to
        ! grid box means.
        do l = 1, land_field
          lake_evap(map_2d(1,ainfo%land_index(l))) = real(fluxes%lake_evap(l), r_def) * &
                                   tile_fraction(map_tile(1,ainfo%land_index(l))+lake-1)
        end do

        ! diagnostics
        if (.not. associated(snomlt_surf_htf, empty_real_data) ) then
          do i = 1, seg_len
            snomlt_surf_htf(map_2d(1,i)) = sf_diag%snomlt_surf_htf(i,1)
          end do
        end if

        if (.not. associated(soil_evap, empty_real_data) ) then
          do i = 1, seg_len
            soil_evap(map_2d(1,i)) = fluxes%esoil_ij_soilt(i,1,1)
          end do
        end if

        if (.not. associated(soil_surf_ht_flux, empty_real_data) ) then
          do i = 1, seg_len
            soil_surf_ht_flux(map_2d(1,i)) = coast%surf_ht_flux_land_ij(i,1)
          end do
        end if

        ! Convert fields to specific quantities, as happens in imp_solver in UM
        if (sf_diag%sq1p5) then
          do i = 1, seg_len
            sf_diag%q1p5m(i,1)=sf_diag%q1p5m(i,1)/ &
                 (1.0_r_def+sf_diag%q1p5m(i,1)+qcf_latest(i,1))
            sf_diag%q1p5m_ssi(i,1)=sf_diag%q1p5m_ssi(i,1)/ &
                 (1.0_r_def+sf_diag%q1p5m_ssi(i,1)+qcf_latest(i,1))
          end do
          do n = 1, ntiles
            do m = 1, ainfo%surft_pts(n)
              l = ainfo%surft_index(m,n)
              i = ainfo%land_index(l)
              sf_diag%q1p5m_surft(l,n) = sf_diag%q1p5m_surft(l,n)/             &
                  (1.0_r_def + sf_diag%q1p5m_surft(l,n)+qcf_latest(i,1) )
            end do
          end do
        end if

        do i = 1, seg_len
          ! Dummy values as unused in ls_cld for lowest level
          cumulus(i,1) = .false.
          ntml(i,1) = 1_i_def
          ! Critical relative humidity
          rhcpt(i,1) = rh_crit_wth(map_wth(1,i) + 1)
        end do

        ! Grid box mean screen level diagnostics
        if (.not. associated(t1p5m, empty_real_data) .or.                      &
             .not. associated(q1p5m, empty_real_data) .or.                     &
             .not. associated(rh1p5m, empty_real_data) .or.                    &
             .not. associated(qcl1p5m, empty_real_data) ) then

          call ls_cld(                                                         &
               forcing%pstar_ij, rhcpt, 1, 1, seg_len, 1, ntml, cumulus,       &
               .false., sf_diag%t1p5m, work_2d_1, sf_diag%q1p5m, qcf_latest,   &
               qcl1p5m_loc, work_2d_2, work_2d_3, error_code )
        end if

        if (.not. associated(t1p5m, empty_real_data) ) then
          do i = 1, seg_len
            t1p5m(map_2d(1,i)) = sf_diag%t1p5m(i,1)
          end do
        end if
        if (.not. associated(q1p5m, empty_real_data) ) then
          do i = 1, seg_len
            q1p5m(map_2d(1,i)) = sf_diag%q1p5m(i,1)
          end do
        end if
        if (.not. associated(qcl1p5m, empty_real_data) ) then
          do i = 1, seg_len
            qcl1p5m(map_2d(1,i)) = qcl1p5m_loc(i,1)
          end do
        end if

        if (.not. associated(rh1p5m, empty_real_data) ) then
          ! qsat needed since q1p5m always a specific humidity
          call qsat(work_2d_1,sf_diag%t1p5m,forcing%pstar_ij,pdims%i_end,pdims%j_end)
          do i = 1, seg_len
            rh1p5m(map_2d(1,i)) = max(0.0_r_def, sf_diag%q1p5m(i,1)) * 100.0_r_def / work_2d_1(i,1)
          end do
        end if

        ! Sea and sea-ice screen level diagnostics
        if (.not. associated(t1p5m_ssi, empty_real_data) .or.                  &
             .not. associated(q1p5m_ssi, empty_real_data) .or.                 &
             .not. associated(rh1p5m_ssi, empty_real_data) .or.                &
             .not. associated(qcl1p5m_ssi, empty_real_data) ) then

          call ls_cld(                                                         &
               forcing%pstar_ij, rhcpt, 1, 1, seg_len, 1, ntml, cumulus,       &
               .false., sf_diag%t1p5m_ssi, work_2d_1, sf_diag%q1p5m_ssi,       &
               qcf_latest, qcl1p5m_loc, work_2d_2, work_2d_3, error_code )
        end if

        if (.not. associated(t1p5m_ssi, empty_real_data) ) then
          do i = 1, seg_len
            t1p5m_ssi(map_2d(1,i)) = sf_diag%t1p5m_ssi(i,1)
          end do
        end if
        if (.not. associated(q1p5m_ssi, empty_real_data) ) then
          do i = 1, seg_len
            q1p5m_ssi(map_2d(1,i)) = sf_diag%q1p5m_ssi(i,1)
          end do
        end if
        if (.not. associated(qcl1p5m_ssi, empty_real_data) ) then
          do i = 1, seg_len
            qcl1p5m_ssi(map_2d(1,i)) = qcl1p5m_loc(i,1)
          end do
        end if

        if (.not. associated(rh1p5m_ssi, empty_real_data) ) then
          ! qsat needed since q1p5m always a specific humidity
          call qsat(work_2d_1,sf_diag%t1p5m_ssi,forcing%pstar_ij,pdims%i_end,pdims%j_end)
          do i = 1, seg_len
            rh1p5m_ssi(map_2d(1,i)) = max(0.0_r_def, sf_diag%q1p5m_ssi(i,1)) * 100.0_r_def / work_2d_1(i,1)
          end do
        end if

        ! Land screen level diagnostics
        if (.not. associated(t1p5m_land, empty_real_data) .or.                 &
             .not. associated(q1p5m_land, empty_real_data) .or.                &
             .not. associated(rh1p5m_land, empty_real_data) .or.               &
             .not. associated(qcl1p5m_land, empty_real_data) ) then

          allocate(t1p5m_land_loc(seg_len,1))
          allocate(q1p5m_land_loc(seg_len,1))
          do i = 1, seg_len
            t1p5m_land_loc(i,1) = 0.0_r_def
            q1p5m_land_loc(i,1) = 0.0_r_def
          end do
          do n = 1, n_land_tile
            do l = 1, land_field
              t1p5m_land_loc(ainfo%land_index(l),1) = &
                   t1p5m_land_loc(ainfo%land_index(l),1) + &
                   ainfo%frac_surft(l,n)*sf_diag%t1p5m_surft(l, n)
              q1p5m_land_loc(ainfo%land_index(l),1) = &
                   q1p5m_land_loc(ainfo%land_index(l),1) + &
                   ainfo%frac_surft(l,n)*sf_diag%q1p5m_surft(l, n)
            end do
          end do

          call ls_cld(                                                         &
               forcing%pstar_ij, rhcpt, 1, 1, seg_len, 1, ntml, cumulus,       &
               .false., t1p5m_land_loc, work_2d_1, q1p5m_land_loc,             &
               qcf_latest, qcl1p5m_loc, work_2d_2, work_2d_3, error_code )
        end if

        if (.not. associated(t1p5m_land, empty_real_data) ) then
          do i = 1, seg_len
            t1p5m_land(map_2d(1,i)) = t1p5m_land_loc(i,1)
          end do
        end if
        if (.not. associated(q1p5m_land, empty_real_data) ) then
          do i = 1, seg_len
            q1p5m_land(map_2d(1,i)) = q1p5m_land_loc(i,1)
          end do
        end if
        if (.not. associated(qcl1p5m_land, empty_real_data) ) then
          do i = 1, seg_len
            qcl1p5m_land(map_2d(1,i)) = qcl1p5m_loc(i,1)
          end do
        end if

        if (.not. associated(rh1p5m_land, empty_real_data) ) then
          ! qsat needed since q1p5m always a specific humidity
          call qsat(work_2d_1,t1p5m_land_loc,forcing%pstar_ij,pdims%i_end,pdims%j_end)
          do i = 1, seg_len
            rh1p5m_land(map_2d(1,i)) = max(0.0_r_def, q1p5m_land_loc(i,1)) * 100.0_r_def / work_2d_1(i,1)
          end do
        end if

        if (allocated(t1p5m_land_loc)) deallocate(t1p5m_land_loc)
        if (allocated(q1p5m_land_loc)) deallocate(q1p5m_land_loc)

        if (.not. associated(t1p5m_surft, empty_real_data) ) then
          do n = 1, n_land_tile
            do l = 1, land_field
              t1p5m_surft(map_tile(1,ainfo%land_index(l))+n-1) = &
                   real(sf_diag%t1p5m_surft(l, n), r_def)
            end do
          end do
        end if

        if (.not. associated(q1p5m_surft, empty_real_data) ) then
          do n = 1, n_land_tile
            do l = 1, land_field
              q1p5m_surft(map_tile(1,ainfo%land_index(l))+n-1) = &
                   real(sf_diag%q1p5m_surft(l, n), r_def)
            end do
          end do
        end if

        if (.not. associated(latent_heat, empty_real_data) ) then
          do n = 1, n_land_tile
            do l = 1, land_field
              latent_heat(map_tile(1,ainfo%land_index(l))+n-1) = &
                   real(fluxes%le_surft(l, n), r_def)
            end do
          end do
          do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
            do i = 1, seg_len
              if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_def) then
                latent_heat(map_tile(1,i)+n-1) = (lc + lf) *                  &
                     tile_moisture_flux(map_tile(1,i)+n-1)
              end if
            end do
          end do
        end if

        if (.not. associated(surf_sw_net, empty_real_data) ) then
          do n = 1, n_land_tile
            do l = 1, land_field
              surf_sw_net(map_tile(1,ainfo%land_index(l))+n-1) =             &
                   real(sw_down_surf(map_2d(1,ainfo%land_index(l))), r_um) - &
                   real(sw_up_tile(map_tile(1,ainfo%land_index(l))+n-1), r_um)
            end do
          end do
          do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
            do i = 1, seg_len
              if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_def) then
                surf_sw_net(map_tile(1,i)+n-1) =                             &
                   real(sw_down_surf(map_2d(1,i)), r_um) -                   &
                   real(sw_up_tile(map_tile(1,i)+n-1), r_um)
              end if
            end do
          end do
          do i = 1, seg_len
            if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then
              surf_sw_net(map_tile(1,i)+first_sea_tile-1) =                    &
                   real(sw_down_surf(map_2d(1,i)), r_um) -                     &
                   real(sw_up_tile(map_tile(1,i)+first_sea_tile-1), r_um)
            end if
          end do
        end if

        if (.not. associated(surf_radnet, empty_real_data) ) then
            do n = 1, n_land_tile
              do l = 1, land_field
                surf_radnet(map_tile(1,ainfo%land_index(l))+n-1) = &
                     real(fluxes%radnet_surft(l, n), r_def)
              end do
            end do
          do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
            i_sice = n - first_sea_ice_tile + 1
            do i = 1, seg_len
              if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_def) then
                surf_radnet(map_tile(1,i)+n-1) =                              &
                     radnet_sice(i,1,i_sice) /                                &
                     ainfo%ice_fract_ncat_sicat(i,1,i_sice)
              end if
            end do
          end do
          do i = 1, seg_len
            if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then
              surf_radnet(map_tile(1,i)+first_sea_tile-1) =                    &
                   real(sw_down_surf(map_2d(1,i)), r_um) -                     &
                   real(sw_up_tile(map_tile(1,i)+first_sea_tile-1), r_um) +    &
                   emis_sea * (forcing%lw_down_ij(i,1) - sbcon * coast%tstar_sea_ij(i,1) ** 4.0_r_def)
            end if
          end do
        end if

        if (.not. associated(surf_lw_up, empty_real_data) ) then
          do n = 1, n_land_tile
            do l = 1, land_field
              surf_lw_up(map_tile(1,ainfo%land_index(l))+n-1) = &
                   real(sf_diag%lw_up_surft(l, n), r_def)
            end do
          end do
          do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
            i_sice = n - first_sea_ice_tile + 1
            do i = 1, seg_len
              if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_def) then
                surf_lw_up(map_tile(1,i)+n-1) =                                &
                     sf_diag%lw_up_sice_weighted_cat(i,1,i_sice) /             &
                     ainfo%ice_fract_ncat_sicat(i,1,i_sice)
              end if
            end do
          end do
          do i = 1, seg_len
            if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then
              surf_lw_up(map_tile(1,i)+first_sea_tile-1) =                     &
                   (1.0_r_def - emis_sea) * forcing%lw_down_ij(i,1) +          &
                   emis_sea * sbcon * coast%tstar_sea_ij(i,1) ** 4.0_r_def
            end if
          end do
        end if

        if (.not. associated(surf_lw_down, empty_real_data) ) then
          do n = 1, n_land_tile
            do l = 1, land_field
              surf_lw_down(map_tile(1,ainfo%land_index(l))+n-1) = &
                   real(sf_diag%lw_down_surft(l, n), r_def)
            end do
          end do
          do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
            do i = 1, seg_len
              if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_def) then
                surf_lw_down(map_tile(1,i)+n-1) = forcing%lw_down_ij(i,1)
              end if
            end do
          end do
          do i = 1, seg_len
            if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then
              surf_lw_down(map_tile(1,i)+first_sea_tile-1) = forcing%lw_down_ij(i,1)
            end if
          end do
        end if

      else !loop = 2

        if (flux_bc_opt == specified_fluxes_only) then
          ! tstar on sea surfaces will not have been updated so, to try and maintain
          ! consistency between tstar and the supplied surface flux, update estimate
          ! of tstar with the increment to level 1 temperature
          ! N.B. dtl1_2d is now the final level 1 T increment, not the first
          ! guess as it was on loop 1
          do i = 1, seg_len
            tile_temperature(map_tile(1,i)+first_sea_tile-1) =                 &
                 tile_temperature(map_tile(1,i)+first_sea_tile-1)+dtl1_2d(map_2d(1,i))
          end do
        end if

        ! Sea-ice bulk temperature
        do n = 1, n_sea_ice_tile
          do i = 1, seg_len
            sea_ice_temperature(map_sice(1,i)+n-1) = real(progs%ti_sicat(i, 1, n), r_def)
          end do
        end do

        ! Sea tile fluxes
        do i = 1, seg_len
          if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then
            tile_heat_flux(map_tile(1,i)+first_sea_tile-1) = ( ftl(i,1) -      &
                 tile_heat_flux(map_tile(1,i)+first_sea_tile-1) )              &
                 / tile_fraction(map_tile(1,i)+first_sea_tile-1)
            tile_moisture_flux(map_tile(1,i)+first_sea_tile-1) = ( fqw(i,1) -  &
                 tile_moisture_flux(map_tile(1,i)+first_sea_tile-1) )          &
                 / tile_fraction(map_tile(1,i)+first_sea_tile-1)
          else
            tile_heat_flux(map_tile(1,i)+first_sea_tile-1) = 0.0_r_def
            tile_moisture_flux(map_tile(1,i)+first_sea_tile-1) = 0.0_r_def
          end if
          snowice_melt(map_tile(1,i)+first_sea_tile-1) = 0.0_r_def
        end do

        if (.not. associated(latent_heat, empty_real_data) ) then
          do i = 1, seg_len
            if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then
              latent_heat(map_tile(1,i)+first_sea_tile-1) = lc *              &
                   tile_moisture_flux(map_tile(1,i)+first_sea_tile-1)
            end if
          end do
        end if

      end if!loop=2

    endif  ! outer = outer_iterations

    ! deallocate diagnostics deallocated in atmos_physics2
    call dealloc_sf_imp(sf_diag)
    deallocate(epot_surft)
    deallocate(tscrndcl_surft)
    deallocate(alpha1)
    deallocate(ashtf_prime_surft)
    deallocate(dtstar_surft)
    deallocate(fracaero_t)
    deallocate(fracaero_s)
    deallocate(rhokh_surft)
    deallocate(chr1p5m)
    deallocate(resfs)
    deallocate(canhc_surft)
    deallocate(resft)
    deallocate(flake)
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

    call aero_nullify(aerotype)
    call aero_dealloc(aero_data)

    call jules_vars_dealloc(jules_vars_data)
    call jules_vars_nullify(jules_vars)

    call prognostics_nullify(progs)
    call prognostics_dealloc(progs_data)

    call fluxes_nullify(fluxes)
    call fluxes_dealloc(fluxes_data)

    call wtrac_jls_nullify(wtrac_jls)
    call wtrac_jls_dealloc(wtrac_jls_data)

  end subroutine jules_imp_code

end module jules_imp_kernel_mod
