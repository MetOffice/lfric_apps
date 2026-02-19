!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the explicit UM boundary layer scheme.
module bl_exp_kernel_mod

  use argument_mod,           only : arg_type,                   &
                                     GH_FIELD, GH_REAL,          &
                                     GH_INTEGER,                 &
                                     GH_READ, GH_WRITE,          &
                                     GH_READWRITE, DOMAIN,       &
                                     ANY_DISCONTINUOUS_SPACE_1,  &
                                     ANY_DISCONTINUOUS_SPACE_2,  &
                                     ANY_DISCONTINUOUS_SPACE_3,  &
                                     ANY_DISCONTINUOUS_SPACE_4,  &
                                     ANY_DISCONTINUOUS_SPACE_5
  use constants_mod,          only : i_def, i_um, r_def, r_um, r_bl
  use empty_data_mod,         only : empty_real_data
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type
  use blayer_config_mod,      only : bl_mix_w
  use cloud_config_mod,       only : rh_crit_opt, rh_crit_opt_tke, scheme,    &
                                     scheme_bimodal, scheme_pc2,              &
                                     pc2ini, pc2ini_bimodal,                  &
                                     i_bm_ez_opt, i_bm_ez_opt_entpar
  use microphysics_config_mod, only: turb_gen_mixph, prog_tnuc
  use mixing_config_mod,       only: smagorinsky
  use jules_surface_config_mod, only : formdrag, formdrag_dist_drag

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_exp_kernel_type
    private
    type(arg_type) :: meta_args(93) = (/                                       &
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! theta_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! rho_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! rho_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! u_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! v_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! w_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! velocity_w2v
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_v_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_cl_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_ci_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! dz_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! rdz_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! dtrdz_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! shear
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! delta
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! max_diff_smag
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_fraction
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sd_orog_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! peak_to_trough_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! silhouette_area_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! rhostar
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! recip_l_mo_sea
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! t1_sd
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! q1_sd
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! dtl_mphys
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! dmt_mphys
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! sw_heating_rate
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! lw_heating_rate
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_liquid
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! rh_crit
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! tnuc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINuOUS_SPACE_1),&! tnuc_nlcl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dsldzm
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! mix_len_bm
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! wvar
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! visc_m_blend
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! visc_h_blend
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dw_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! rhokm_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_3),&! surf_interp
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, W3),                       &! rhokh_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! tke_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! ngstress_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! bq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! bt_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, W3),                       &! moist_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, W3),                       &! heat_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! fd_taux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! fd_tauy
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sea_u_current
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sea_v_current
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! lmix_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! gradrinr
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! z0m_eff
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ustar
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! zhsc_2d
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! z_lcl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! inv_depth
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qcl_at_inv_top
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_1),&! shallow_flag
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! uw0_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! vw0_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lcl_height
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! parcel_top
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_1),&! level_parcel_top
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! wstar_2d
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! thv_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! parcel_buoyancy
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qsat_at_lcl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! bl_weight_1dbl
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_4),&! bl_type_ind
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),&! level_ent
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),&! level_ent_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! ent_we_lim
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! ent_t_frac
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! ent_zrzi
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! ent_we_lim_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! ent_t_frac_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! ent_zrzi_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! diag__zht
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1) &! diag__oblen
         /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: bl_exp_code
  end type

  public :: bl_exp_code

contains

  !> @brief Interface to the UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     theta_in_wth           Potential temperature field
  !> @param[in]     rho_in_w3              Density field in density space
  !> @param[in]     rho_in_wth             Density field in theta space
  !> @param[in]     wetrho_in_wth          Wet density field in wth space
  !> @param[in]     exner_in_w3            Exner pressure field in density space
  !> @param[in]     exner_in_wth           Exner pressure field in wth space
  !> @param[in]     u_in_w3                'Zonal' wind in density space
  !> @param[in]     v_in_w3                'Meridional' wind in density space
  !> @param[in]     w_in_wth               'Vertical' wind in theta space
  !> @param[in]     velocity_w2v           Velocity normal to cell top
  !> @param[in]     m_v_n                  Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n                 Cloud liquid mixing ratio at time level n
  !> @param[in]     m_ci_n                 Cloud ice mixing ratio at time level n
  !> @param[in]     height_w3              Height of density space above surface
  !> @param[in]     height_wth             Height of theta space above surface
  !> @param[in]     dz_wth                 Layer depths at wtheta points
  !> @param[in]     rdz_w3                 Inverse Layer depths at w3 points
  !> @param[in]     dtrdz_wth              dt/(rho*r*r*dz) in wth
  !> @param[in]     shear                  3D wind shear on wtheta points
  !> @param[in]     delta                  Edge length on wtheta points
  !> @param[in]     max_diff_smag          Maximum diffusion coefficient allowed
  !> @param[in,out] zh_2d                  Boundary layer depth
  !> @param[in,out] ntml_2d                Number of turbulently mixed levels
  !> @param[in,out] cumulus_2d             Cumulus flag (true/false)
  !> @param[in]     tile_fraction          Surface tile fractions
  !> @param[in]     sd_orog_2d             Standard deviation of orography
  !> @param[in]     peak_to_trough_orog    Half of peak-to-trough height over root(2) of orography
  !> @param[in]     silhouette_area_orog   Silhouette area of orography
  !> @param[in]     tile_temperature       Surface tile temperatures
  !> @param[in]     rhostar_2d             Surface density
  !> @param[in]     recip_l_mo_sea_2d      Inverse Obukhov length over sea only
  !> @param[in]     t1_sd_2d               StDev of level 1 temperature
  !> @param[in]     q1_sd_2d               StDev of level 1 humidity
  !> @param[in]     dtl_mphys              Microphysics liq temperature increment
  !> @param[in]     dmt_mphys              Microphysics total water increment
  !> @param[in]     sw_heating_rate        Shortwave radiation heating rate
  !> @param[in]     lw_heating_rate        Longwave radiation heating rate
  !> @param[in]     cf_bulk                Bulk cloud fraction
  !> @param[in]     cf_liquid              Liquid cloud fraction
  !> @param[in,out] rh_crit                Critical rel humidity
  !> @param[in]     tnuc                   Temperature of nucleation (K)
  !> @param[in,out] tnuc_nlcl              Temperature of nucleation (K) (2D)
  !> @param[in,out] dsldzm                 Liquid potential temperature gradient in wth
  !> @param[in,out] mix_len_bm             Turb length-scale for bimodal in wth
  !> @param[in,out] wvar                   Vertical velocity variance in wth
  !> @param[in,out] visc_m_blend           Blended BL-Smag diffusion coefficient for momentum
  !> @param[in,out] visc_h_blend           Blended BL-Smag diffusion coefficient for scalars
  !> @param[in,out] dw_bl                  Vertical wind increment from BL scheme
  !> @param[in,out] rhokm_bl               Momentum eddy diffusivity on BL levels
  !> @param[in,out] surf_interp            Surface variables for regridding
  !> @param[in,out] rhokh_bl               Heat eddy diffusivity on BL levels
  !> @param[in,out] tke_bl                 Turbulent kinetic energy (m2 s-2)
  !> @param[in,out] ngstress_bl            Non-gradient stress function on BL levels
  !> @param[in,out] bq_bl                  Buoyancy parameter for moisture
  !> @param[in,out] bt_bl                  Buoyancy parameter for heat
  !> @param[in,out] moist_flux_bl          Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl           Vertical heat flux on BL levels
  !> @param[in,out] dtrdz_tq_bl            dt/(rho*r*r*dz) in wth
  !> @param[in,out] fd_taux                'Zonal' momentum stress from form drag
  !> @param[in,out] fd_tauy                'Meridional' momentum stress from form drag
  !> @param[in]     sea_u_current          Ocean surface U current
  !> @param[in]     sea_v_current          Ocean surface V current
  !> @param[in,out] lmix_bl                Turbulence mixing length in wth
  !> @param[in,out] gradrinr               Gradient Richardson number in wth
  !> @param[in]     z0m_eff                Grid mean effective roughness length
  !> @param[in]     ustar                  Friction velocity
  !> @param[in,out] zh_nonloc              Depth of non-local BL scheme
  !> @param[in,out] zhsc_2d                Height of decoupled layer top
  !> @param[in,out] z_lcl                  Height of the LCL (wtheta levels)
  !> @param[in,out] inv_depth              Depth of BL top inversion layer
  !> @param[in,out] qcl_at_inv_top         Cloud water at top of inversion
  !> @param[in,out] shallow_flag           Indicator of shallow convection
  !> @param[in,out] uw0_flux               'Zonal' surface momentum flux
  !> @param[in,out] vw0 flux               'Meridional' surface momentum flux
  !> @param[in,out] lcl_height             Height of lifting condensation level (w3 levels)
  !> @param[in,out] parcel_top             Height of surface based parcel ascent
  !> @param[in,out] level_parcel_top       Model level of parcel_top
  !> @param[in,out] wstar_2d               BL velocity scale
  !> @param[in,out] thv_flux               Surface flux of theta_v
  !> @param[in,out] parcel_buoyancy        Integral of parcel buoyancy
  !> @param[in,out] qsat_at_lcl            Saturation specific hum at LCL
  !> @param[in,out] bl_weight_1dbl         Blending weight to 1D BL scheme in the BL
  !> @param[in,out] bl_type_ind            Diagnosed BL types
  !> @param[in,out] level_ent              Level of surface mixed layer inversion
  !> @param[in,out] level_ent_dsc          Level of decoupled stratocumulus inversion
  !> @param[in,out] ent_we_lim             Rho * entrainment rate at surface ML inversion (kg m-2 s-1)
  !> @param[in,out] ent_t_frac             Fraction of time surface ML inversion is above level
  !> @param[in,out] ent_zrzi               Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in,out] ent_we_lim_dsc         Rho * entrainment rate at DSC inversion (kg m-2 s-1)
  !> @param[in,out] ent_t_frac_dsc         Fraction of time DSC inversion is above level
  !> @param[in,out] ent_zrzi_dsc           Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in,out] zht                    Diagnostic: turb mixing height
  !> @param[in,out] oblen                  Diagnostic: Obukhov length
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
  !> @param[in]     ndf_surf               Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf              Number of unique DOFs for surface variables
  !> @param[in]     map_surf               Dofmap for the cell at the base of the column for surface variables
  !> @param[in]     ndf_bl                 Number of DOFs per cell for BL types
  !> @param[in]     undf_bl                Number of total DOFs for BL types
  !> @param[in]     map_bl                 Dofmap for cell for BL types
  !> @param[in]     ndf_ent                Number of DOFs per cell for entrainment levels
  !> @param[in]     undf_ent               Number of total DOFs for entrainment levels
  !> @param[in]     map_ent                Dofmap for cell for entrainment levels
  subroutine bl_exp_code(nlayers, seg_len,                      &
                         theta_in_wth,                          &
                         rho_in_w3,                             &
                         rho_in_wth,                            &
                         wetrho_in_wth,                         &
                         exner_in_w3,                           &
                         exner_in_wth,                          &
                         u_in_w3,                               &
                         v_in_w3,                               &
                         w_in_wth,                              &
                         velocity_w2v,                          &
                         m_v_n,                                 &
                         m_cl_n,                                &
                         m_ci_n,                                &
                         height_w3,                             &
                         height_wth,                            &
                         dz_wth,                                &
                         rdz_w3,                                &
                         dtrdz_wth,                             &
                         shear,                                 &
                         delta,                                 &
                         max_diff_smag,                         &
                         zh_2d,                                 &
                         ntml_2d,                               &
                         cumulus_2d,                            &
                         tile_fraction,                         &
                         sd_orog_2d,                            &
                         peak_to_trough_orog,                   &
                         silhouette_area_orog,                  &
                         tile_temperature,                      &
                         rhostar_2d,                            &
                         recip_l_mo_sea_2d,                     &
                         t1_sd_2d,                              &
                         q1_sd_2d,                              &
                         dtl_mphys,                             &
                         dmt_mphys,                             &
                         sw_heating_rate,                       &
                         lw_heating_rate,                       &
                         cf_bulk,                               &
                         cf_liquid,                             &
                         rh_crit,                               &
                         tnuc,                                  &
                         tnuc_nlcl,                             &
                         dsldzm,                                &
                         mix_len_bm,                            &
                         wvar,                                  &
                         visc_m_blend,                          &
                         visc_h_blend,                          &
                         dw_bl,                                 &
                         rhokm_bl,                              &
                         surf_interp,                           &
                         rhokh_bl,                              &
                         tke_bl,                                &
                         ngstress_bl,                           &
                         bq_bl,                                 &
                         bt_bl,                                 &
                         moist_flux_bl,                         &
                         heat_flux_bl,                          &
                         dtrdz_tq_bl,                           &
                         fd_taux,                               &
                         fd_tauy,                               &
                         sea_u_current,                         &
                         sea_v_current,                         &
                         lmix_bl,                               &
                         gradrinr,                              &
                         z0m_eff,                               &
                         ustar,                                 &
                         zh_nonloc,                             &
                         zhsc_2d,                               &
                         z_lcl,                                 &
                         inv_depth,                             &
                         qcl_at_inv_top,                        &
                         shallow_flag,                          &
                         uw0_flux,                              &
                         vw0_flux,                              &
                         lcl_height,                            &
                         parcel_top,                            &
                         level_parcel_top,                      &
                         wstar_2d,                              &
                         thv_flux,                              &
                         parcel_buoyancy,                       &
                         qsat_at_lcl,                           &
                         bl_weight_1dbl,                        &
                         bl_type_ind,                           &
                         level_ent,                             &
                         level_ent_dsc,                         &
                         ent_we_lim,                            &
                         ent_t_frac,                            &
                         ent_zrzi,                              &
                         ent_we_lim_dsc,                        &
                         ent_t_frac_dsc,                        &
                         ent_zrzi_dsc,                          &
                         zht,                                   &
                         oblen,                                 &
                         ndf_wth, undf_wth, map_wth,            &
                         ndf_w3, undf_w3, map_w3,               &
                         ndf_2d, undf_2d, map_2d,               &
                         ndf_tile, undf_tile, map_tile,         &
                         ndf_surf, undf_surf, map_surf,         &
                         ndf_bl, undf_bl, map_bl,               &
                         ndf_ent, undf_ent, map_ent)

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use jules_control_init_mod, only: n_surf_tile

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use atm_fields_bounds_mod, only: pdims
    use bl_option_mod, only: alpha_cd, l_noice_in_turb, l_use_surf_in_ri
    use cv_run_mod, only: i_convection_vn, i_convection_vn_6a,               &
                          cldbase_opt_dp, cldbase_opt_md
    use nlsizes_namelist_mod, only: bl_levels
    use planet_constants_mod, only: p_zero, kappa, planet_radius, &
                                    lcrcp => lcrcp_bl, lsrcp => lsrcp_bl
    use timestep_mod, only: timestep

    use free_tracers_inputs_mod,    only: n_wtrac
    use wtrac_atm_step_mod,         only: atm_step_wtrac_type
    use wtrac_bl_mod,               only: bl_wtrac_type

    ! subroutines used
    use atmos_physics2_save_restore_mod, only: ap2_init_conv_diag
    use bl_diags_mod, only: BL_diag, dealloc_bl_imp, dealloc_bl_expl, &
                            alloc_bl_expl
    use conv_diag_6a_mod, only: conv_diag_6a
    use buoy_tq_mod, only: buoy_tq
    use bdy_expl2_mod, only: bdy_expl2
    use tr_mix_mod, only: tr_mix

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth, seg_len)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3, seg_len)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d, seg_len)

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile, seg_len)

    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf, ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: map_surf(ndf_surf, seg_len)
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl, seg_len)
    integer(kind=i_def), intent(in) :: ndf_ent, undf_ent
    integer(kind=i_def), intent(in) :: map_ent(ndf_ent, seg_len)

    real(kind=r_def), dimension(undf_wth), intent(inout):: rh_crit,            &
                                                           dsldzm,             &
                                                           mix_len_bm,         &
                                                           wvar, dw_bl,        &
                                                           visc_h_blend,       &
                                                           visc_m_blend,       &
                                                           rhokm_bl,           &
                                                           tke_bl,             &
                                                           ngstress_bl,        &
                                                           bq_bl, bt_bl,       &
                                                           dtrdz_tq_bl,        &
                                                           lmix_bl,            &
                                                           gradrinr
    real(kind=r_def), dimension(undf_w3),  intent(inout):: rhokh_bl,           &
                                                           moist_flux_bl,      &
                                                           heat_flux_bl,       &
                                                           fd_taux, fd_tauy
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: rho_in_w3,          &
                                                           exner_in_w3,        &
                                                           u_in_w3, v_in_w3,   &
                                                           height_w3, rdz_w3
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           rho_in_wth,         &
                                                           wetrho_in_wth,      &
                                                           exner_in_wth,       &
                                                           w_in_wth,           &
                                                           velocity_w2v,       &
                                                           m_v_n, m_cl_n,      &
                                                           m_ci_n,             &
                                                           height_wth,         &
                                                           dz_wth,             &
                                                           dtrdz_wth,          &
                                                           shear, delta,       &
                                                           max_diff_smag,      &
                                                           dtl_mphys,dmt_mphys,&
                                                           sw_heating_rate,    &
                                                           lw_heating_rate,    &
                                                           cf_bulk, cf_liquid, &
                                                           tnuc
    real(kind=r_def), dimension(undf_2d), intent(inout) :: zh_2d,              &
                                                           zhsc_2d,            &
                                                           z0m_eff,            &
                                                           ustar,              &
                                                           zh_nonloc,          &
                                                           z_lcl,              &
                                                           inv_depth,          &
                                                           qcl_at_inv_top,     &
                                                           uw0_flux,           &
                                                           vw0_flux,           &
                                                           lcl_height,         &
                                                           parcel_top,         &
                                                           wstar_2d,           &
                                                           thv_flux,           &
                                                           parcel_buoyancy,    &
                                                           qsat_at_lcl,        &
                                                           bl_weight_1dbl,     &
                                                           tnuc_nlcl
    integer(kind=i_def), dimension(undf_2d), intent(inout) :: ntml_2d,         &
                                                              cumulus_2d,      &
                                                              shallow_flag,    &
                                                              level_parcel_top
    real(kind=r_def), dimension(undf_2d), intent(in)    :: recip_l_mo_sea_2d,  &
                                                           rhostar_2d,         &
                                                           t1_sd_2d, q1_sd_2d, &
                                                           sea_u_current,      &
                                                           sea_v_current

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)
    real(kind=r_def), intent(in) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(in) :: sd_orog_2d(undf_2d)
    real(kind=r_def), intent(in) :: peak_to_trough_orog(undf_2d)
    real(kind=r_def), intent(in) :: silhouette_area_orog(undf_2d)

    integer(kind=i_def), dimension(undf_bl), intent(inout) :: bl_type_ind
    real(kind=r_def), dimension(undf_surf), intent(inout)  :: surf_interp
    integer(kind=i_def), dimension(undf_2d), intent(inout) :: level_ent
    integer(kind=i_def), dimension(undf_2d), intent(inout) :: level_ent_dsc
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_we_lim
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_t_frac
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_zrzi
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_we_lim_dsc
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_t_frac_dsc
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_zrzi_dsc

    real(kind=r_def), pointer, intent(inout) :: zht(:)
    real(kind=r_def), pointer, intent(inout) :: oblen(:)
    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    integer(i_def) :: k, i, l, n, land_field

    ! local switches and scalars
    integer(i_um) :: error_code
    real(r_bl) :: weight1, weight2, weight3
    logical :: l_spec_z0, l_cape_opt
    logical, parameter :: l_extra_call = .false.

    ! profile fields from level 1 upwards
    real(r_bl), dimension(seg_len,1,nlayers) :: rho_dry, z_rho, z_theta,     &
         bulk_cloud_fraction, rho_wet_tq, u_p, v_p, rhcpt, theta,            &
         p_rho_levels, exner_rho_levels, tgrad_bm, mix_len_tmp,              &
         exner_theta_levels,                                                 &
         bulk_cf_conv, qcf_conv, r_rho_levels, visc_h, visc_m, rneutml_sq,    &
         tnuc_new

    ! profile field on boundary layer levels
    real(r_bl), dimension(seg_len,1,bl_levels) :: fqw, ftl, rhokh, bq_gb,    &
         bt_gb, dtrdz_charney_grid, rdz_charney_grid, rhokm_mix,             &
         temperature, rho_mix_tq, dzl_charney, qw, tl, bt, bq,               &
         bt_cld, bq_cld, a_qs, a_dqsdt, dqsdt, rhokm, tau_fd_x, tau_fd_y, rdz

    real(r_um), dimension(seg_len,1,bl_levels) :: w_mixed, w_flux

    ! profile fields from level 2 upwards
    real(r_bl), dimension(seg_len,1,2:nlayers+1) :: bl_w_var

    real(r_bl), dimension(seg_len,1,2:bl_levels) :: f_ngstress

    ! profile fields from level 0 upwards
    real(r_bl), dimension(seg_len,1,0:nlayers) :: p_theta_levels, etadot, w, &
         q, qcl, qcf, r_theta_levels

    ! profile fields with a hard-wired 2
    real(r_bl), dimension(seg_len,1,2,bl_levels) :: rad_hr, micro_tends

    ! single level real fields
    real(r_bl), dimension(seg_len,1) :: p_star, tstar, zh_prev, zlcl, zhpar, &
         zh, dzh, wstar, wthvs, u_0_p, v_0_p, zlcl_uv, qsat_lcl, delthvu,    &
         bl_type_1, bl_type_2, bl_type_3, bl_type_4, bl_type_5, bl_type_6,   &
         bl_type_7, uw0, vw0, zhnl, rhostar,                                 &
         recip_l_mo_sea, flandg, t1_sd, q1_sd, qcl_inv_top,                  &
         fb_surf, rib_gb, z0m_eff_gb, zhsc, ustargbm, cos_theta_latitude,    &
         max_diff, delta_smag, tnuc_nlcl_um
    real(r_um), dimension(seg_len,1) :: surf_dep_flux, zeroes

    real(r_bl), dimension(seg_len,1,3) :: t_frac, t_frac_dsc, we_lim, &
         we_lim_dsc, zrzi, zrzi_dsc

    ! single level integer fields
    integer(i_um), dimension(seg_len,1) :: ntml, ntpar, kent, kent_dsc

    ! single level logical fields
    logical, dimension(seg_len,1) :: land_sea_mask, cumulus, l_shallow

    ! fields on land points
    real(r_bl), dimension(:), allocatable :: sil_orog_land_gb, ho2r2_orog_gb, &
         sd_orog

    ! integer fields on land points
    integer, dimension(:), allocatable :: land_index

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    integer(i_um), parameter :: nscmdpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs)=.false.

    real(r_bl), dimension(seg_len,1,nlayers) :: rho_wet

    real(r_bl), dimension(seg_len,1,0:nlayers) :: conv_prog_precip

    real(r_bl), dimension(seg_len,1) :: z0h_scm, z0m_scm, w_max, ql_ad,      &
         cin_undilute, cape_undilute, entrain_coef, ustar_in, g_ccp, h_ccp,  &
         ccp_strength, cu_over_orog, shallowc, flux_e, flux_h,               &
         z0msea, tstar_sea, tstar_land, ice_fract, tstar_sice

    integer(i_um), dimension(seg_len,1) :: nlcl, conv_type, nbdsc, ntdsc

    logical, dimension(seg_len,1) :: no_cumulus, l_congestus, l_congestus2,  &
         l_mid

    ! Water tracer fields which are not currently used but are required
    ! UM routine
    type (atm_step_wtrac_type), dimension(n_wtrac) :: wtrac_as
    type (bl_wtrac_type), dimension(n_wtrac) :: wtrac_bl

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    !-----------------------------------------------------------------------
    ! diagnostic flags
    error_code=0

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------

    ! Land fraction
    land_field = 0
    do i = 1, seg_len
      flandg(i,1) = surf_interp(map_surf(1,i)+0)
      if (flandg(i,1) > 0.0_r_bl) then
        land_field = land_field + 1
      end if
      fb_surf(i,1) = surf_interp(map_surf(1,i)+6)
    end do

    allocate(land_index(land_field))
    l = 0
    do i = 1, seg_len
      if (flandg(i,1) > 0.0_r_bl) then
        l = l+1
        land_index(l) = i
      end if
    end do

    if (l_use_surf_in_ri) then
      do i = 1, seg_len
        tstar(i,1) = 0.0_r_bl
        do n = 1, n_surf_tile
          if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_bl) then
            tstar(i,1) = tstar(i,1) + tile_temperature(map_tile(1,i)+n-1) * tile_fraction(map_tile(1,i)+n-1)
          end if
        end do
      end do
    end if

    allocate(sd_orog(land_field))
    allocate(ho2r2_orog_gb(land_field))
    allocate(sil_orog_land_gb(land_field))

    do l = 1, land_field
      ! Standard deviation of orography
      sd_orog(l) = real(sd_orog_2d(map_2d(1,land_index(l))), r_bl)
      ! Half of peak-to-trough height over root(2) of orography (ho2r2_orog_gb)
      ho2r2_orog_gb(l) = real(peak_to_trough_orog(map_2d(1,land_index(l))), r_bl)
      sil_orog_land_gb(l) = real(silhouette_area_orog(map_2d(1,land_index(l))), r_bl)
    end do

    ! Information passed from Jules explicit
    do i = 1, seg_len
      ustargbm(i,1) = ustar(map_2d(1,i))
      rhostar(i,1) = rhostar_2d(map_2d(1,i))
      recip_l_mo_sea(i,1) = recip_l_mo_sea_2d(map_2d(1,i))
      rib_gb(i,1) = gradrinr(map_wth(1,i))
      z0m_eff_gb(i,1) = z0m_eff(map_2d(1,i))
      ftl(i,1,1) = heat_flux_bl(map_w3(1,i))
      fqw(i,1,1) = moist_flux_bl(map_w3(1,i))
      rhokh(i,1,1) = rhokh_bl(map_w3(1,i))
      rhokm(i,1,1) = rhokm_bl(map_wth(1,i))
      t1_sd(i,1) = t1_sd_2d(map_2d(1,i))
      q1_sd(i,1) = q1_sd_2d(map_2d(1,i))
    end do

    if (prog_tnuc) then
      ! Use tnuc from LFRic and map onto tnuc_new for UM to be passed to conv_diag_6a
      do k = 1, nlayers
        do i = 1, seg_len
          tnuc_new(i,1,k) = real(tnuc(map_wth(1,i) + k),kind=r_bl)
        end do ! i
      end do ! k
    end if

    !-----------------------------------------------------------------------
    ! assuming map_wth(1,i) points to level 0
    ! and map_w3(1,i) points to level 1
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      do k = 0, nlayers
        ! w wind on theta levels
        w(i,1,k) = w_in_wth(map_wth(1,i) + k)
        ! height of theta levels from centre of planet
        r_theta_levels(i,1,k) = height_wth(map_wth(1,i) + k) + planet_radius
        p_theta_levels(i,1,k) = p_zero*(exner_in_wth(map_wth(1,i) + k))**(1.0_r_def/kappa)
      end do
      do k = 1, nlayers
        exner_theta_levels(i,1,k) = exner_in_wth(map_wth(1,i) + k)
        ! potential temperature on theta levels
        theta(i,1,k) = theta_in_wth(map_wth(1,i) + k)
        ! wet density on theta and rho levels
        rho_wet_tq(i,1,k) = wetrho_in_wth(map_wth(1,i) + k)
        ! dry density on rho levels
        rho_dry(i,1,k) = rho_in_w3(map_w3(1,i) + k-1)
        ! pressure on rho levels
        p_rho_levels(i,1,k) = p_zero*(exner_in_w3(map_w3(1,i) + k-1))**(1.0_r_def/kappa)
        ! exner pressure on rho levels
        exner_rho_levels(i,1,k) = exner_in_w3(map_w3(1,i) + k-1)
        ! u wind on rho levels
        u_p(i,1,k) = u_in_w3(map_w3(1,i) + k-1)
        ! v wind on rho levels
        v_p(i,1,k) = v_in_w3(map_w3(1,i) + k-1)
        ! height of rho levels from centre of planet
        r_rho_levels(i,1,k) = height_w3(map_w3(1,i) + k-1) + planet_radius
        ! height of levels above surface
        z_rho(i,1,k) = r_rho_levels(i,1,k)-r_theta_levels(i,1,0)
        z_theta(i,1,k) = r_theta_levels(i,1,k)-r_theta_levels(i,1,0)
        ! water vapour mixing ratio
        q(i,1,k) = m_v_n(map_wth(1,i) + k)
        ! cloud liquid mixing ratio
        qcl(i,1,k) = m_cl_n(map_wth(1,i) + k)
        ! cloud ice mixing ratio
        qcf_conv(i,1,k) = m_ci_n(map_wth(1,i) + k)
        bulk_cf_conv(i,1,k) = cf_bulk(map_wth(1,i) + k)
      end do
    end do
    if (l_noice_in_turb) then
      do i = 1, seg_len
        do k = 1, nlayers
          qcf(i,1,k) = 0.0_r_bl
          bulk_cloud_fraction(i,1,k) = cf_liquid(map_wth(1,i) + k)
        end do
      end do
    else
      do i = 1, seg_len
        do k = 1, nlayers
          qcf(i,1,k) = m_ci_n(map_wth(1,i) + k)
          bulk_cloud_fraction(i,1,k) = cf_bulk(map_wth(1,i) + k)
        end do
      end do
    end if

    if ( smagorinsky ) then
      do i = 1, seg_len
        delta_smag(i,1) = delta(map_wth(1,i))
        max_diff(i,1) = max_diff_smag(map_wth(1,i))
        do k = 1, nlayers
          visc_m(i,1,k) = shear(map_wth(1,i) + k)
          visc_h(i,1,k) = shear(map_wth(1,i) + k)
        end do
      end do
    end if

    ! Set this to 1 to account for quasi-uniform grid
    do i = 1, seg_len
      cos_theta_latitude(i,1) = 1.0_r_um
    end do

    do i = 1, seg_len
      ! surface pressure
      p_star(i,1) = p_theta_levels(i,1,0)
      do k = 1, nlayers
        ! computational vertical velocity
        etadot(i,1,k) = velocity_w2v(map_wth(1,i) + k) / z_theta(i,1,nlayers)
      end do
      ! surface currents
      u_0_p(i,1) = 0.0_r_bl
      v_0_p(i,1) = 0.0_r_bl
    end do

    do i = 1, seg_len
      ! surface currents
      u_0_p(i,1) = sea_u_current(map_2d(1,i))
      v_0_p(i,1) = sea_v_current(map_2d(1,i))
    end do

    !-----------------------------------------------------------------------
    ! Things saved from one timestep to the next
    !-----------------------------------------------------------------------
    ! previous BL height
    do i = 1, seg_len
      zh(i,1) = zh_2d(map_2d(1,i))
      zh_prev(i,1) = zh(i,1)
    end do

    !-----------------------------------------------------------------------
    ! Things saved from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      do k = 1, bl_levels
        ! microphysics tendancy terms
        micro_tends(i,1,1,k) = dtl_mphys(map_wth(1,i)+k)/timestep
        micro_tends(i,1,2,k) = dmt_mphys(map_wth(1,i)+k)/timestep
        ! radiation tendancy terms
        rad_hr(i,1,1,k) = lw_heating_rate(map_wth(1,i)+k)
        rad_hr(i,1,2,k) = sw_heating_rate(map_wth(1,i)+k)
        ! temperature
        temperature(i,1,k) = theta(i,1,k) * exner_theta_levels(i,1,k)
        tl(i,1,k) = temperature(i,1,k) - lcrcp*qcl(i,1,k) - lsrcp*qcf(i,1,k)
        qw(i,1,k) = q(i,1,k) + qcl(i,1,k) + qcf(i,1,k)
      end do
    end do

    !-----------------------------------------------------------------------
    ! Boundary layer diagnostics
    !-----------------------------------------------------------------------

    ! needed to ensure zht is saved if wanted
    bl_diag%l_zht     = .not. associated(zht, empty_real_data)
    bl_diag%l_oblen   = .not. associated(oblen, empty_real_data)
    bl_diag%l_weight1d = .true.
    call alloc_bl_expl(bl_diag, .true.)

    bl_diag%l_tke      = .true.
    bl_diag%l_elm3d    = .true.
    bl_diag%l_gradrich = .true.
    allocate(BL_diag%tke(seg_len,1,bl_levels))
    allocate(BL_diag%elm3d(seg_len,1,bl_levels))
    allocate(BL_diag%gradrich(seg_len,1,bl_levels))
    do k = 1, bl_levels
      do i = 1, seg_len
        bl_diag%tke(i,1,k) = 0.0_r_bl
        bl_diag%elm3d(i,1,k) = 0.0_r_bl
        bl_diag%gradrich(i,1,k) = 0.0_r_bl
      end do
    end do

    ! Calculate vertical differences
    do i = 1, seg_len
      dzl_charney(i,1,1) = 2.0_r_bl * (r_theta_levels(i,1,1) - r_theta_levels(i,1,0))
      do k = 2, bl_levels
        dzl_charney(i,1,k) = dz_wth(map_wth(1,i) + k)
        rdz(i,1,k) = 1.0_r_bl/dz_wth(map_wth(1,i) + k-1)
      end do
      do k = 1, bl_levels
        rdz_charney_grid(i,1,k) = rdz_w3(map_w3(1,i) + k-1)
        rho_mix_tq(i,1,k) = rho_in_wth(map_wth(1,i) + k)
        dtrdz_charney_grid(i,1,k) = dtrdz_wth(map_wth(1,i) + k) / rho_mix_tq(i,1,k)
      end do
    end do

    call buoy_tq (                                                             &
      ! IN dimensions/logicals
      bl_levels,                                                               &
      ! IN fields
      p_theta_levels,temperature,q,qcf,qcl,bulk_cloud_fraction,                &
      ! OUT fields
      bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                       &
      )

    ! Use  convection switches to decide the value of  L_cape_opt
    if (i_convection_vn == i_convection_vn_6a ) then
      L_cape_opt = ( (cldbase_opt_dp == 3) .or. (cldbase_opt_md == 3) .or. &
                     (cldbase_opt_dp == 4) .or. (cldbase_opt_md == 4) .or. &
                     (cldbase_opt_dp == 5) .or. (cldbase_opt_md == 5) .or. &
                     (cldbase_opt_dp == 6) .or. (cldbase_opt_md == 6) )
    else
      L_cape_opt = .false.
    end if

    call ap2_init_conv_diag( 1, seg_len, ntml, ntpar, nlcl, cumulus,        &
        l_shallow, l_mid, delthvu, ql_ad, zhpar, dzh, qcl_inv_top,          &
        zlcl, zlcl_uv, conv_type, no_cumulus, w_max, w, L_cape_opt)

    do i = 1, seg_len
      qsat_lcl(i,1) = 0.0_r_bl
    end do
    call conv_diag_6a(                                                  &
    !     IN Parallel variables
            seg_len, 1                                                  &
    !     IN model dimensions.
          , bl_levels, p_rho_levels, p_theta_levels(1,1,1)              &
          , exner_rho_levels, rho_wet, rho_wet_tq, z_theta, z_rho       &
          , r_theta_levels, r_rho_levels                                &
    !     IN Model switches
          , l_extra_call, no_cumulus                                    &
    !     IN cloud data
          , qcf_conv, qcl(1,1,1), bulk_cf_conv                          &
    !     IN everything not covered so far :
          , p_star, q(1,1,1), theta, exner_theta_levels, u_p, v_p       &
          , u_0_p, v_0_p, tstar_land, tstar_sea, tstar_sice, z0msea     &
          , flux_e, flux_h, ustar_in, L_spec_z0, z0m_scm, z0h_scm       &
          , tstar, land_sea_mask, flandg, ice_fract, w, w_max           &
          , conv_prog_precip, g_ccp, h_ccp, ccp_strength                &
    !     IN surface fluxes
          , fb_surf, ustargbm                                           &
    !     SCM Diagnostics (dummy values in full UM)
          , nSCMDpkgs,L_SCMDiags                                        &
    !     OUT data required elsewhere in UM system :
          , zh,zhpar,dzh,qcl_inv_top,zlcl,zlcl_uv,delthvu,ql_ad, ntml   &
          , ntpar,nlcl, cumulus,l_shallow,l_congestus,l_congestus2      &
          , conv_type, CIN_undilute,CAPE_undilute, wstar, wthvs         &
          , entrain_coef, qsat_lcl, Error_code, tnuc_new, tnuc_nlcl_um )

    call bdy_expl2 (                                                           &
    ! IN values defining vertical grid of model atmosphere :
      bl_levels,p_theta_levels,land_field,land_index,                          &
      r_theta_levels, r_rho_levels, cos_theta_latitude,                        &
    ! IN U, V and W momentum fields.
      u_p,v_p, u_0_p, v_0_p,                                                   &
    ! IN from other part of explicit boundary layer code
      rho_dry,rho_wet_tq,rho_mix_tq,dzl_charney,rdz,rdz_charney_grid,          &
      z_theta,z_rho,rhostar,bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,      &
      dqsdt,recip_l_mo_sea, flandg, rib_gb, sil_orog_land_gb,z0m_eff_gb,       &
    ! IN cloud/moisture data :
      bulk_cloud_fraction,q,qcf,qcl,temperature,qw,tl,                         &
    ! IN everything not covered so far :
      rad_hr,micro_tends,fb_surf,ustargbm,p_star,tstar,                        &
      zh_prev, zhpar,zlcl,ho2r2_orog_gb,sd_orog,wtrac_as,                      &
    ! 2 IN 3 INOUT for Smagorinsky
      delta_smag, max_diff, rneutml_sq, visc_m, visc_h,                        &
    ! SCM Diagnostics (dummy values in full UM) & stash diag
      nSCMDpkgs,L_SCMDiags,BL_diag,                                            &
    ! INOUT variables
      zh,dzh,ntml,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,w,etadot,        &
      t1_sd,q1_sd,wtrac_bl,                                                    &
    ! OUT New variables for message passing
      tau_fd_x, tau_fd_y, f_ngstress,                                          &
    ! OUT Diagnostic not requiring STASH flags :
      zhnl,shallowc,cu_over_orog,                                              &
      bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6, bl_type_7,  &
    ! OUT data for turbulent generation of mixed-phase cloud:
      bl_w_var,                                                                &
    ! OUT data required for tracer mixing :
      kent, we_lim, t_frac, zrzi, kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,  &
    ! OUT data required elsewhere in UM system :
      zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0,rhcpt, tgrad_bm, mix_len_tmp        &
      )

    if (prog_tnuc) then
      ! Use tnuc_nlcl_um from conv_diag_6a (UM) and map onto tnuc_nlcl for
      ! LFRic to then be passed out to
      do i = 1, seg_len
        tnuc_nlcl(map_2d(1,i)) = real(tnuc_nlcl_um(i,1),kind=r_def)
      end do
    end if

    if (bl_mix_w) then

      ! Interpolate rhokm from theta-levels to rho-levels, as is done for
      ! rhokh in bdy_expl2.
      ! NOTE: rhokm is defined on theta-levels with the k-indexing offset by
      ! 1 compared to the rest of the UM (k=1 is the surface).

      ! Bottom model-level is surface in both arrays, so no interp needed
      ! (for rhokh_mix, this is done in the JULES routine sf_impl2_jls).
      do i = 1, seg_len
        rhokm_mix(i,1,1) = rhokm(i,1,1)
      end do
      do k = 2, bl_levels-1
        do i = 1, seg_len
          weight1 = r_theta_levels(i,1,k) - r_theta_levels(i,1,k-1)
          weight2 = r_theta_levels(i,1,k) - r_rho_levels(i,1,k)
          weight3 = r_rho_levels(i,1,k)   - r_theta_levels(i,1,k-1)
          rhokm_mix(i,1,k) = (weight3/weight1) * rhokm(i,1,k+1) &
                           + (weight2/weight1) * rhokm(i,1,k)
          ! Scale exchange coefficients by 1/dz factor, as is done for
          ! rhokh_mix in bdy_impl4
          ! (note this doesn't need to be done for the surface exchange coef)
          rhokm_mix(i,1,k) = rhokm_mix(i,1,k) * rdz_charney_grid(i,1,k)
        end do
      end do
      k = bl_levels
      do i = 1, seg_len
        weight1 = r_theta_levels(i,1,k) - r_theta_levels(i,1,k-1)
        weight2 = r_theta_levels(i,1,k) - r_rho_levels(i,1,k)
        ! Assume rhokm(BL_LEVELS+1) is zero
        rhokm_mix(i,1,k) = (weight2/weight1) * rhokm(i,1,k)
        ! Scale exchange coefficients by 1/dz factor, as is done for
        ! rhokh_mix in bdy_impl4
        rhokm_mix(i,1,k) = rhokm_mix(i,1,k) * rdz_charney_grid(i,1,k)
        zeroes(i,1) = 0.0_r_um
      end do

      do k = 1, bl_levels
        do i = 1, seg_len
          w_mixed(i,1,k) = w(i,1,k)
        end do
      end do

      call  tr_mix (                                                           &
           ! IN fields
           real(r_theta_levels,r_um), real(r_rho_levels,r_um), pdims,          &
           bl_levels, alpha_cd,                                                &
           real(rhokm_mix(1:seg_len,1:1,2:bl_levels),r_um),                    &
           real(rhokm_mix(1:seg_len,1:1,1),r_um),                              &
           real(dtrdz_charney_grid,r_um), zeroes, zeroes, kent,                &
           real(we_lim,r_um), real(t_frac,r_um), real(zrzi,r_um), kent_dsc,    &
           real(we_lim_dsc,r_um), real(t_frac_dsc,r_um), real(zrzi_dsc,r_um),  &
           real(zhnl,r_um), real(zhsc,r_um), real(z_rho,r_um),                 &
           ! INOUT / OUT fields
           w_mixed, w_flux, surf_dep_flux                                      &
           )

      do k = 1, bl_levels
        do i = 1, seg_len
          dw_bl(map_wth(1,i)+k) = w_mixed(i,1,k) - w(i,1,k)
        end do
      end do

    end if

    ! 2D variables that need interpolating to cell faces
    do i = 1, seg_len
      surf_interp(map_surf(1,i)+5) = zhnl(i,1)
    end do

    do k=1,bl_levels
      do i = 1, seg_len
        bq_bl(map_wth(1,i) + k-1) = bq_gb(i,1,k)
        bt_bl(map_wth(1,i) + k-1) = bt_gb(i,1,k)
        dtrdz_tq_bl(map_wth(1,i) + k) = dtrdz_charney_grid(i,1,k)
      end do
    end do

    do i = 1, seg_len
      level_ent(map_2d(1,i)) = int( kent(i,1), i_def )
      level_ent_dsc(map_2d(1,i)) = int( kent_dsc(i,1), i_def )
    end do
    do k = 1, 3
      do i = 1, seg_len
        ent_we_lim(map_ent(1,i) + k - 1) = real( we_lim(i,1,k), r_def )
        ent_t_frac(map_ent(1,i) + k - 1) = real( t_frac(i,1,k), r_def )
        ent_zrzi(map_ent(1,i) + k - 1) = real( zrzi(i,1,k), r_def )
        ent_we_lim_dsc(map_ent(1,i) + k - 1) = real( we_lim_dsc(i,1,k), r_def )
        ent_t_frac_dsc(map_ent(1,i) + k - 1) = real( t_frac_dsc(i,1,k), r_def )
        ent_zrzi_dsc(map_ent(1,i) + k - 1) = real( zrzi_dsc(i,1,k), r_def )
      end do
    end do

    if (formdrag == formdrag_dist_drag) then
      do k = 1, bl_levels
        do i = 1, seg_len
          ! These fields will be passed to set wind, which maps w3 (cell centre)
          ! to w2 (cell face) vectors. However, they are actually defined in
          ! wtheta (cell top centre) and need mapping to fd1 (cell top edge).
          ! Set wind will therefore work correctly, but the indexing is shifted
          ! by half a level in the vertical for the input & output
          fd_taux(map_w3(1,i) + k-1) = tau_fd_x(i,1,k)
          fd_tauy(map_w3(1,i) + k-1) = tau_fd_y(i,1,k)
        end do
      end do
      do i = 1, seg_len
        do k=bl_levels+1,nlayers
          fd_taux(map_w3(1,i) + k-1) = 0.0_r_def
          fd_tauy(map_w3(1,i) + k-1) = 0.0_r_def
        end do
      end do
    end if

    do i = 1, seg_len
      gradrinr(map_wth(1,i)) = rib_gb(i,1)
    end do
    do k = 2, bl_levels
      do i = 1, seg_len
        rhokm_bl(map_wth(1,i) + k-1) = rhokm(i,1,k)
        rhokh_bl(map_w3(1,i) + k-1) = rhokh(i,1,k)
        moist_flux_bl(map_w3(1,i) + k-1) = fqw(i,1,k)
        heat_flux_bl(map_w3(1,i) + k-1) = ftl(i,1,k)
        gradrinr(map_wth(1,i) + k-1) = BL_diag%gradrich(i,1,k)
        lmix_bl(map_wth(1,i) + k-1)  = BL_diag%elm3d(i,1,k)
        ngstress_bl(map_wth(1,i) + k-1) = f_ngstress(i,1,k)
        tke_bl(map_wth(1,i) + k-1) = BL_diag%tke(i,1,k)
      end do
    end do

    do i = 1, seg_len
      bl_weight_1dbl(map_2d(1,i)) = BL_diag%weight1d(i,1,2)
      zh_nonloc(map_2d(1,i)) = zhnl(i,1)
      z_lcl(map_2d(1,i)) = real(zlcl(i,1), r_def)
      inv_depth(map_2d(1,i)) = real(dzh(i,1), r_def)
      qcl_at_inv_top(map_2d(1,i)) = real(qcl_inv_top(i,1), r_def)
      if ( l_shallow(i,1) ) then
        shallow_flag(map_2d(1,i)) = 1_i_def
      else
        shallow_flag(map_2d(1,i)) = 0_i_def
      end if
      uw0_flux(map_2d(1,i)) = uw0(i,1)
      vw0_flux(map_2d(1,i)) = vw0(i,1)
      lcl_height(map_2d(1,i)) = zlcl_uv(i,1)
      parcel_top(map_2d(1,i)) = zhpar(i,1)
      level_parcel_top(map_2d(1,i)) = ntpar(i,1)
      wstar_2d(map_2d(1,i)) = wstar(i,1)
      thv_flux(map_2d(1,i)) = wthvs(i,1)
      parcel_buoyancy(map_2d(1,i)) = delthvu(i,1)
      qsat_at_lcl(map_2d(1,i)) = qsat_lcl(i,1)

      bl_type_ind(map_bl(1,i)+0) = bl_type_1(i,1)
      bl_type_ind(map_bl(1,i)+1) = bl_type_2(i,1)
      bl_type_ind(map_bl(1,i)+2) = bl_type_3(i,1)
      bl_type_ind(map_bl(1,i)+3) = bl_type_4(i,1)
      bl_type_ind(map_bl(1,i)+4) = bl_type_5(i,1)
      bl_type_ind(map_bl(1,i)+5) = bl_type_6(i,1)
      bl_type_ind(map_bl(1,i)+6) = bl_type_7(i,1)
    end do

    ! update blended Smagorinsky diffusion coefficients only if using Smagorinsky scheme
    if ( smagorinsky ) then
      do i = 1, seg_len
        visc_m_blend(map_wth(1,i)) = visc_m(i,1,1)
        visc_h_blend(map_wth(1,i)) = visc_h(i,1,1)
      end do
      do k = 1, bl_levels-1
        do i = 1, seg_len
          visc_m_blend(map_wth(1,i) + k) = visc_m(i,1,k)
          visc_h_blend(map_wth(1,i) + k) = visc_h(i,1,k)
        end do
      end do
      do i = 1, seg_len
        visc_m_blend(map_wth(1,i) + bl_levels) = visc_m(i,1,bl_levels-1)
        visc_h_blend(map_wth(1,i) + bl_levels) = visc_h(i,1,bl_levels-1)
      end do
      do k = bl_levels+1, nlayers
        do i = 1, seg_len
          visc_m_blend(map_wth(1,i) + k) = 0.0_r_def
          visc_h_blend(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    endif

    ! update BL prognostics
    do i = 1, seg_len
      zh_2d(map_2d(1,i))     = zh(i,1)
      zhsc_2d(map_2d(1,i))   = zhsc(i,1)
      ntml_2d(map_2d(1,i))   = ntml(i,1)
      if (cumulus(i,1)) then
        cumulus_2d(map_2d(1,i)) = 1_i_def
      else
        cumulus_2d(map_2d(1,i)) = 0_i_def
      endif
    end do

    if (rh_crit_opt == rh_crit_opt_tke) then
      do i = 1, seg_len
        rh_crit(map_wth(1,i)) = real(rhcpt(i,1,1), r_def)
      end do
      do k = 1, nlayers
        do i = 1, seg_len
          rh_crit(map_wth(1,i)+k) = real(rhcpt(i,1,k), r_def)
        end do
      end do
    end if

    ! Liquid temperature gradient and turbulent length-scale
    ! for bimodal cloud scheme
    if (scheme == scheme_bimodal .or. &
         (scheme == scheme_pc2 .and. pc2ini == pc2ini_bimodal ) ) then
      if (i_bm_ez_opt == i_bm_ez_opt_entpar) then
        ! Length-scale used for entraining parcel mode construction method
        do k = 1, nlayers
          do i = 1, seg_len
            mix_len_bm(map_wth(1,i)+k) = mix_len_tmp(i,1,k)
          end do
        end do
      else
        ! SL-gradient used for stable-layer mode construction method
        do k = 1, nlayers
          do i = 1, seg_len
            dsldzm(map_wth(1,i)+k) = tgrad_bm(i,1,k)
          end do
        end do
      end if
    end if
    if (scheme == scheme_bimodal .or. turb_gen_mixph .or. &
         (scheme == scheme_pc2 .and. pc2ini == pc2ini_bimodal ) ) then
      do k = 2, nlayers+1
        do i = 1, seg_len
          wvar(map_wth(1,i)+k-1) = bl_w_var(i,1,k)
        end do
      end do
    end if

    if (.not. associated(zht, empty_real_data) ) then
      do i = 1, seg_len
        zht(map_2d(1,i)) = BL_diag%zht(i,1)
      end do
    end if
    if (.not. associated(oblen, empty_real_data) ) then
      do i = 1, seg_len
        oblen(map_2d(1,i)) = BL_diag%oblen(i,1)
      end do
    end if

    ! deallocate diagnostics deallocated in atmos_physics2
    call dealloc_bl_expl(bl_diag)
    deallocate(BL_diag%gradrich)
    deallocate(BL_diag%elm3d)
    deallocate(BL_diag%tke)
    deallocate(land_index)
    deallocate(sd_orog)
    deallocate(ho2r2_orog_gb)
    deallocate(sil_orog_land_gb)

  end subroutine bl_exp_code

end module bl_exp_kernel_mod
