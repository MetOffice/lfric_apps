!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM High Order Turbulence Closure Scheme.
module bl_hoc_exp_kernel_mod

  use argument_mod,              only: arg_type,                   &
                                       GH_FIELD, GH_REAL,          &
                                       GH_INTEGER,                 &
                                       GH_READ, GH_WRITE,          &
                                       GH_READWRITE, DOMAIN,       &
                                       ANY_DISCONTINUOUS_SPACE_1,  &
                                       ANY_DISCONTINUOUS_SPACE_2,  &
                                       ANY_DISCONTINUOUS_SPACE_3,  &
                                       ANY_DISCONTINUOUS_SPACE_4
  use constants_mod,             only: i_def, i_um, r_def, r_um, r_bl
  use empty_data_mod,            only: empty_real_data
  use fs_continuity_mod,         only: W3, Wtheta
  use kernel_mod,                only: kernel_type
  use mixing_config_mod,         only: smagorinsky, fullstress, leonard_tke
  use blayer_config_mod,         only: shcu_buoy, bdy_tke, bdy_tke_deardorff
  use mym_option_mod,            only: tke_levels
  use microphysics_config_mod,   only: prog_tnuc
  use jules_surface_config_mod,  only: formdrag, formdrag_dist_drag

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_hoc_exp_kernel_type
    private
    type(arg_type) :: meta_args(84) = (/                                       &
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! theta_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! rho_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! rho_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! u_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! v_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! w_in_wth
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
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_fraction
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sd_orog_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! peak_to_trough_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! silhouette_area_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_liquid
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! tnuc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINuOUS_SPACE_1),&! tnuc_nlcl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! visc_m_blend
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! visc_h_blend
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dw_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! rhokm_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_3),&! surf_interp
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, W3),                       &! rhokh_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! tke_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! tsq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! qsq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! cov_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! zhpar_shcu_2d
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! rhogamu_w3
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! rhogamv_w3
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! leonard_klm_tke
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! leonard_klh_tke
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! bq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! bt_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, W3),                       &! moist_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, W3),                       &! heat_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! fd_taux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! fd_tauy
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sea_u_current
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sea_v_current
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! master_length
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! gradrinr
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! rhogamu_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! rhogamv_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! rhogamt_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! rhogamq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! tke_shr_prod
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! tke_boy_prod
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! tke_dissp
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! sm25
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! sh25
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dbdz
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dvdzm
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! z0m_eff
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ustar
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
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_4) &! bl_type_ind
         /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: bl_hoc_exp_code
  end type

  public :: bl_hoc_exp_code

contains

  !> @brief Interface to the UM High Order Turbulence Closure Scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP25
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
  !> @param[in,out] zh_2d                  Boundary layer depth
  !> @param[in,out] ntml_2d                Number of turbulently mixed levels
  !> @param[in,out] cumulus_2d             Cumulus flag (true/false)
  !> @param[in]     tile_fraction          Surface tile fractions
  !> @param[in]     sd_orog_2d             Standard deviation of orography
  !> @param[in]     peak_to_trough_orog    Half of peak-to-trough height over root(2) of orography
  !> @param[in]     silhouette_area_orog   Silhouette area of orography
  !> @param[in]     tile_temperature       Surface tile temperatures
  !> @param[in]     cf_bulk                Bulk cloud fraction
  !> @param[in]     cf_liquid              Liquid cloud fraction
  !> @param[in]     tnuc                   Temperature of nucleation (K)
  !> @param[in,out] tnuc_nlcl              Temperature of nucleation (K) (2D)
  !> @param[in,out] visc_m_blend           Blended BL-Smag diffusion coefficient for momentum
  !> @param[in,out] visc_h_blend           Blended BL-Smag diffusion coefficient for scalars
  !> @param[in,out] dw_bl                  Vertical wind increment from BL scheme
  !> @param[in,out] rhokm_bl               Momentum eddy diffusivity on BL levels
  !> @param[in,out] surf_interp            Surface variables for regridding
  !> @param[in,out] rhokh_bl               Heat eddy diffusivity on BL levels
  !> @param[in,out] tke_bl                 Turbulent kinetic energy (m2 s-2)
  !> @param[in,out] tsq_bl                 Self covariant of tl'
  !> @param[in,out] qsq_bl                 Self covariant of qw'
  !> @param[in,out] cov_bl                 Correlation between tl' and qw'
  !> @param[in,out] zhpar_shcu_2d          Mixed layer height for non-gradient buoyancy flux
  !> @param[in,out] rhogamu_w3             Counter Gradient Flux Term for U
  !> @param[in,out] rhogamv_w3             Counter Gradient Flux Term for V
  !> @param[in,out] leonard_klm_tke        Leonard term coefficient for momentum
  !> @param[in,out] leonard_klh_tke        Leonard term coefficient for heat
  !> @param[in,out] bq_bl                  Buoyancy parameter for moisture
  !> @param[in,out] bt_bl                  Buoyancy parameter for heat
  !> @param[in,out] moist_flux_bl          Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl           Vertical heat flux on BL levels
  !> @param[in,out] dtrdz_tq_bl            dt/(rho*r*r*dz) in wth
  !> @param[in,out] fd_taux                'Zonal' momentum stress from form drag
  !> @param[in,out] fd_tauy                'Meridional' momentum stress from form drag
  !> @param[in]     sea_u_current          Ocean surface U current
  !> @param[in]     sea_v_current          Ocean surface V current
  !> @param[in,out] master_length          Turbulence length scale in wth
  !> @param[in,out] gradrinr               Gradient Richardson number in wth
  !> @param[in,out] rhogamu_bl             Counter Gradient Flux Term for U
  !> @param[in,out] rhogamv_bl             Counter Gradient Flux Term for V
  !> @param[in,out] rhogamt_bl             Counter Gradient Flux Term for tl
  !> @param[in,out] rhogamq_bl             Counter Gradient Flux Term for qt
  !> @param[in,out] tke_shr_prod           Production rate of TKE by shear
  !> @param[in,out] tke_boy_prod           Production rate of TKE by buoyancy
  !> @param[in,out] tke_dissp              Dissipation rate of TKE
  !> @param[in,out] sm25                   Stability function for momentum
  !> @param[in,out] sh25                   Stability function for scalar
  !> @param[in,out] dbdz                   Vertical gradient of buoyancy
  !> @param[in,out] dvdzm                  Modulus of wind shear
  !> @param[in]     z0m_eff                Grid mean effective roughness length
  !> @param[in]     ustar                  Friction velocity
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
  !> @param[in,out] bl_type_ind            Diagnosed BL types
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
  subroutine bl_hoc_exp_code(nlayers, seg_len,                    &
                             theta_in_wth,                        &
                             rho_in_w3,                           &
                             rho_in_wth,                          &
                             wetrho_in_wth,                       &
                             exner_in_w3,                         &
                             exner_in_wth,                        &
                             u_in_w3,                             &
                             v_in_w3,                             &
                             w_in_wth,                            &
                             m_v_n,                               &
                             m_cl_n,                              &
                             m_ci_n,                              &
                             height_w3,                           &
                             height_wth,                          &
                             dz_wth,                              &
                             rdz_w3,                              &
                             dtrdz_wth,                           &
                             shear,                               &
                             delta,                               &
                             zh_2d,                               &
                             ntml_2d,                             &
                             cumulus_2d,                          &
                             tile_fraction,                       &
                             sd_orog_2d,                          &
                             peak_to_trough_orog,                 &
                             silhouette_area_orog,                &
                             tile_temperature,                    &
                             cf_bulk,                             &
                             cf_liquid,                           &
                             tnuc,                                &
                             tnuc_nlcl,                           &
                             visc_m_blend,                        &
                             visc_h_blend,                        &
                             dw_bl,                               &
                             rhokm_bl,                            &
                             surf_interp,                         &
                             rhokh_bl,                            &
                             tke_bl,                              &
                             tsq_bl,                              &
                             qsq_bl,                              &
                             cov_bl,                              &
                             zhpar_shcu_2d,                       &
                             rhogamu_w3,                          &
                             rhogamv_w3,                          &
                             leonard_klm_tke,                     &
                             leonard_klh_tke,                     &
                             bq_bl,                               &
                             bt_bl,                               &
                             moist_flux_bl,                       &
                             heat_flux_bl,                        &
                             dtrdz_tq_bl,                         &
                             fd_taux,                             &
                             fd_tauy,                             &
                             sea_u_current,                       &
                             sea_v_current,                       &
                             master_length,                       &
                             gradrinr,                            &
                             rhogamu_bl,                          &
                             rhogamv_bl,                          &
                             rhogamt_bl,                          &
                             rhogamq_bl,                          &
                             tke_shr_prod,                        &
                             tke_boy_prod,                        &
                             tke_dissp,                           &
                             sm25,                                &
                             sh25,                                &
                             dbdz,                                &
                             dvdzm,                               &
                             z0m_eff,                             &
                             ustar,                               &
                             z_lcl,                               &
                             inv_depth,                           &
                             qcl_at_inv_top,                      &
                             shallow_flag,                        &
                             uw0_flux,                            &
                             vw0_flux,                            &
                             lcl_height,                          &
                             parcel_top,                          &
                             level_parcel_top,                    &
                             wstar_2d,                            &
                             thv_flux,                            &
                             parcel_buoyancy,                     &
                             qsat_at_lcl,                         &
                             bl_type_ind,                         &
                             ndf_wth, undf_wth, map_wth,          &
                             ndf_w3, undf_w3, map_w3,             &
                             ndf_2d, undf_2d, map_2d,             &
                             ndf_tile, undf_tile, map_tile,       &
                             ndf_surf, undf_surf, map_surf,       &
                             ndf_bl, undf_bl, map_bl)

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

    use wtrac_atm_step_mod,         only: atm_step_wtrac_type

    ! subroutines used
    use atmos_physics2_save_restore_mod, only: ap2_init_conv_diag
    use bl_diags_mod, only: BL_diag, dealloc_bl_imp, dealloc_bl_expl, &
                            alloc_bl_expl
    use conv_diag_6a_mod, only: conv_diag_6a
    use buoy_tq_mod, only: buoy_tq
    use bdy_expl2_1a_mod, only: bdy_expl2_1a
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

    real(kind=r_def), dimension(undf_wth), intent(inout):: dw_bl,              &
                                                           visc_h_blend,       &
                                                           visc_m_blend,       &
                                                           rhokm_bl,           &
                                                           tke_bl,             &
                                                           tsq_bl,             &
                                                           qsq_bl,             &
                                                           cov_bl,             &
                                                           leonard_klm_tke,    &
                                                           bq_bl, bt_bl,       &
                                                           dtrdz_tq_bl,        &
                                                           gradrinr
    real(kind=r_def), dimension(undf_w3),  intent(inout):: rho_in_w3,          &
                                                           rhokh_bl,           &
                                                           moist_flux_bl,      &
                                                           heat_flux_bl,       &
                                                           fd_taux, fd_tauy,   &
                                                           rhogamu_w3,         &
                                                           rhogamv_w3,         &
                                                           leonard_klh_tke
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: exner_in_w3,        &
                                                           u_in_w3, v_in_w3,   &
                                                           height_w3, rdz_w3
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           rho_in_wth,         &
                                                           wetrho_in_wth,      &
                                                           exner_in_wth,       &
                                                           w_in_wth,           &
                                                           m_v_n, m_cl_n,      &
                                                           m_ci_n,             &
                                                           height_wth,         &
                                                           dz_wth,             &
                                                           dtrdz_wth,          &
                                                           shear, delta,       &
                                                           cf_bulk, cf_liquid, &
                                                           tnuc
    real(kind=r_def), dimension(undf_2d), intent(inout) :: zh_2d,              &
                                                           zhpar_shcu_2d,      &
                                                           z0m_eff,            &
                                                           ustar,              &
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
                                                           tnuc_nlcl
    integer(kind=i_def), dimension(undf_2d), intent(inout) :: ntml_2d,         &
                                                              cumulus_2d,      &
                                                              shallow_flag,    &
                                                              level_parcel_top
    real(kind=r_def), dimension(undf_2d), intent(in)    :: sea_u_current,      &
                                                           sea_v_current

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)
    real(kind=r_def), intent(in) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(in) :: sd_orog_2d(undf_2d)
    real(kind=r_def), intent(in) :: peak_to_trough_orog(undf_2d)
    real(kind=r_def), intent(in) :: silhouette_area_orog(undf_2d)

    integer(kind=i_def), dimension(undf_bl), intent(inout) :: bl_type_ind
    real(kind=r_def), dimension(undf_surf), intent(inout)  :: surf_interp

    real(kind=r_def), pointer, intent(inout) :: master_length(:)
    real(kind=r_def), pointer, intent(inout) :: rhogamu_bl(:)
    real(kind=r_def), pointer, intent(inout) :: rhogamv_bl(:)
    real(kind=r_def), pointer, intent(inout) :: rhogamt_bl(:)
    real(kind=r_def), pointer, intent(inout) :: rhogamq_bl(:)
    real(kind=r_def), pointer, intent(inout) :: tke_shr_prod(:)
    real(kind=r_def), pointer, intent(inout) :: tke_boy_prod(:)
    real(kind=r_def), pointer, intent(inout) :: tke_dissp(:)
    real(kind=r_def), pointer, intent(inout) :: sm25(:)
    real(kind=r_def), pointer, intent(inout) :: sh25(:)
    real(kind=r_def), pointer, intent(inout) :: dbdz(:)
    real(kind=r_def), pointer, intent(inout) :: dvdzm(:)
    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    integer(i_def) :: k, i, l, n, land_field

    ! local switches and scalars
    integer(i_um) :: error_code
    real(r_bl) :: weight1, weight2, weight3
    logical :: l_spec_z0, l_cape_opt

    ! profile fields from level 1 upwards
    real(r_bl), dimension(seg_len,1,nlayers) :: z_rho, z_theta,              &
         bulk_cloud_fraction, rho_wet_tq, u_p, v_p, theta,                   &
         p_rho_levels, exner_rho_levels,                                     &
         exner_theta_levels,                                                 &
         bulk_cf_conv, qcf_conv, tnuc_new
    ! Single precision is not accurate enough for distance from centre of planet
    real(r_um), dimension(seg_len,1,nlayers) :: r_rho_levels

    ! profile field on boundary layer levels
    real(r_bl), dimension(seg_len,1,bl_levels) :: fqw, ftl, rhokh, bq_gb,    &
         bt_gb, rdz_charney_grid,                                            &
         temperature, qw, tl, bt, bq,                                        &
         bt_cld, bq_cld, a_qs, a_dqsdt, dqsdt, rhokm, tau_fd_x, tau_fd_y, rdz, &
         shear_3d, visc_m, visc_h, tke_trb, tsq_trb, qsq_trb, cov_trb

    real(r_bl), dimension(seg_len,1,bl_levels+1) :: rho_mix

    real(r_um), dimension(seg_len,1,bl_levels) :: w_mixed, w_flux,           &
                                                  rhokm_mix,                 &
                                                  dtrdz_charney_grid

    real(r_bl), dimension(seg_len,1,2:bl_levels) :: rhogamu, rhogamv

    real(r_bl), dimension(seg_len,1,bl_levels,2) :: leonard_kl_tke

    ! profile fields from level 0 upwards
    real(r_bl), dimension(seg_len,1,0:nlayers) :: p_theta_levels, w,         &
         q, qcl, qcf
    ! Single precision is not accurate enough for distance from centre of planet
    real(r_um), dimension(seg_len,1,0:nlayers) :: r_theta_levels

    ! single level real fields
    real(r_bl), dimension(seg_len,1) :: p_star, tstar, zh_prev, zlcl, zhpar, &
         zh, dzh, wstar, wthvs, u_0_p, v_0_p, zlcl_uv, qsat_lcl, delthvu,    &
         bl_type_1, bl_type_2, bl_type_3, bl_type_4, bl_type_5, bl_type_6,   &
         bl_type_7, uw0, vw0,                                                &
         h_blend_orog, flandg, qcl_inv_top,                                  &
         fb_surf, rib_gb, z0m_eff_gb, zhsc, ustargbm,                        &
         delta_smag, tnuc_nlcl_um
    real(r_um), dimension(seg_len,1) :: surf_dep_flux

    real(r_bl), dimension(seg_len,1,3) :: t_frac, t_frac_dsc, we_lim,        &
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
         z0msea, tstar_sea, tstar_land, ice_fract, tstar_sice,               &
         zhpar_shcu

    ! single level integer fields
    integer(i_um), dimension(seg_len,1) :: nlcl, conv_type, nbdsc, ntdsc

    ! single level logical fields
    logical, dimension(seg_len,1) :: no_cumulus, l_congestus, l_congestus2,  &
         l_mid

    integer, dimension(seg_len,1) :: kent_dummy
    real(r_um), dimension(seg_len,1) :: zeroes_2d
    real(r_um), dimension(seg_len,1,3) :: zeroes_ent
    logical, parameter :: l_extra_call = .false.

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
      tstar = 0.0_r_bl
      do i = 1, seg_len
        do n = 1, n_surf_tile
          if (tile_fraction(map_tile(1,i)+n-1) > 0.0_r_bl) then
            tstar(i,1) = tstar(i,1) + tile_temperature(map_tile(1,i)+n-1) *    &
                         tile_fraction(map_tile(1,i)+n-1)
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
      rib_gb(i,1) = gradrinr(map_wth(1,i))
      z0m_eff_gb(i,1) = z0m_eff(map_2d(1,i))
      ftl(i,1,1) = heat_flux_bl(map_w3(1,i))
      fqw(i,1,1) = moist_flux_bl(map_w3(1,i))
      rhokh(i,1,1) = rhokh_bl(map_w3(1,i))
      rhokm(i,1,1) = rhokm_bl(map_wth(1,i))
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
        z_rho(i,1,k) = height_w3(map_w3(1,i) + k-1) - height_wth(map_wth(1,i))
        z_theta(i,1,k) = height_wth(map_wth(1,i) + k) - height_wth(map_wth(1,i))
        ! water vapour mixing ratio
        q(i,1,k) = m_v_n(map_wth(1,i) + k)
        ! cloud liquid mixing ratio
        qcl(i,1,k) = m_cl_n(map_wth(1,i) + k)
        ! cloud ice mixing ratio
        qcf_conv(i,1,k) = m_ci_n(map_wth(1,i) + k)
        bulk_cf_conv(i,1,k) = cf_bulk(map_wth(1,i) + k)
        if (l_noice_in_turb) then
          qcf(i,1,k) = 0.0_r_bl
          bulk_cloud_fraction(i,1,k) = cf_liquid(map_wth(1,i) + k)
        else
          qcf(i,1,k) = m_ci_n(map_wth(1,i) + k)
          bulk_cloud_fraction(i,1,k) = cf_bulk(map_wth(1,i) + k)
        end if
      end do

      do k = 1, bl_levels+1
        rho_mix(i,1,k) = rho_in_w3(map_w3(1,i) + k-1)
      end do

      do k = 1, bl_levels
        temperature(i,1,k) = theta_in_wth(map_wth(1,i) + k) * &
                             exner_in_wth(map_wth(1,i) + k)
        tl(i,1,k) = temperature(i,1,k) - lcrcp*qcl(i,1,k) - lsrcp*qcf(i,1,k)
        qw(i,1,k) = q(i,1,k) + qcl(i,1,k) + qcf(i,1,k)
        rdz_charney_grid(i,1,k) = rdz_w3(map_w3(1,i) + k-1)
        dtrdz_charney_grid(i,1,k) = dtrdz_wth(map_wth(1,i) + k) /              &
                                      rho_in_wth(map_wth(1,i) + k)
        tke_trb(i,1,k) = tke_bl(map_wth(1,i) + k-1)
      end do

      if (bdy_tke /= bdy_tke_deardorff) then
        do k = 1, bl_levels
          tsq_trb(i,1,k) = tsq_bl(map_wth(1,i) + k-1)
          qsq_trb(i,1,k) = qsq_bl(map_wth(1,i) + k-1)
          cov_trb(i,1,k) = cov_bl(map_wth(1,i) + k-1)
        end do
      end if

      do k = 2, bl_levels
        rdz(i,1,k) = 1.0_r_bl/dz_wth(map_wth(1,i) + k-1)
      end do

      ! surface pressure
      p_star(i,1) = p_theta_levels(i,1,0)
      ! surface currents
      u_0_p(i,1) = sea_u_current(map_2d(1,i))
      v_0_p(i,1) = sea_v_current(map_2d(1,i))
      ! previous BL height
      zh(i,1) = zh_2d(map_2d(1,i))
      zh_prev(i,1) = zh(i,1)
    end do

    if (shcu_buoy) then
      do i = 1, seg_len
        zhpar_shcu(i,1) = zhpar_shcu_2d(map_2d(1,i))
      end do
    end if

    if ( smagorinsky ) then
      do i = 1, seg_len
        delta_smag(i,1) = delta(map_wth(1,i))
        do k = 1, bl_levels
          shear_3d(i,1,k) = shear(map_wth(1,i) + k)
        end do
      end do
    end if

    !-----------------------------------------------------------------------
    ! Boundary layer diagnostics
    !-----------------------------------------------------------------------

    bl_diag%l_rhogamt      = .not. associated(rhogamt_bl, empty_real_data)
    bl_diag%l_rhogamq      = .not. associated(rhogamq_bl, empty_real_data)
    bl_diag%l_elm          = .not. associated(master_length, empty_real_data)
    bl_diag%l_tke_shr_prod = .not. associated(tke_shr_prod, empty_real_data)
    bl_diag%l_tke_boy_prod = .not. associated(tke_boy_prod, empty_real_data)
    bl_diag%l_tke_dissp    = .not. associated(tke_dissp, empty_real_data)
    bl_diag%l_sm           = .not. associated(sm25, empty_real_data)
    bl_diag%l_sh           = .not. associated(sh25, empty_real_data)
    bl_diag%l_dbdz         = .not. associated(dbdz, empty_real_data)
    bl_diag%l_dvdzm        = .not. associated(dvdzm, empty_real_data)

    call alloc_bl_expl(bl_diag, .true.)

    bl_diag%l_gradrich = .true.
    allocate(BL_diag%gradrich(seg_len,1,bl_levels))
    bl_diag%gradrich = 0.0_r_bl

    !----------------------------------------------------------------!
    ! Run boundary layer scheme
    !----------------------------------------------------------------!
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

    qsat_lcl = 0.0_r_bl
    call conv_diag_6a(                                                  &
    !     IN Parallel variables
            seg_len, 1                                                  &
    !     IN model dimensions.
          , bl_levels, p_rho_levels, p_theta_levels(1,1,1)              &
          , exner_rho_levels, rho_wet, rho_wet_tq, z_theta, z_rho       &
          , r_theta_levels                                              &
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

    if (prog_tnuc) then
      ! Use tnuc_nlcl_um from conv_diag_6a (UM) and map onto tnuc_nlcl for
      ! LFRic to then be passed out to
      do i = 1, seg_len
        tnuc_nlcl(map_2d(1,i)) = real(tnuc_nlcl_um(i,1),kind=r_def)
      end do
    end if

    call bdy_expl2_1a (                                                        &
    ! IN values defining vertical grid of model atmosphere :
      bl_levels,p_theta_levels,land_field,land_index,                          &
    ! IN U, V and W momentum fields.
      u_p,v_p, u_0_p, v_0_p,                                                   &
    ! IN variables for TKE scheme
      p_star,p_rho_levels,                                                     &
    ! IN from other part of explicit boundary layer code
      rho_mix,rho_wet_tq,rdz,rdz_charney_grid,                                 &
      z_theta,z_rho,bt,bt_gb,bq_gb,                                            &
      flandg, rib_gb, sil_orog_land_gb,z0m_eff_gb,                             &
    ! IN cloud/moisture data:
      q,qcf,qcl,temperature,qw,tl,                                             &
    ! IN everything not covered so far :
      fb_surf,ustargbm,                                                        &
      zh_prev,ho2r2_orog_gb,sd_orog,                                           &
    ! 2 IN for Smagorinsky
      delta_smag, shear_3d,                                                    &
    ! stash diag
      BL_diag,                                                                 &
    ! INOUT variables
      zh,ntml,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,                     &
      tke_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                          &
    ! OUT New variables for message passing
      tau_fd_x, tau_fd_y, visc_m, visc_h, rhogamu, rhogamv,                    &
    ! OUT Diagnostic not requiring STASH flags :
      shallowc,cu_over_orog,                                                   &
      bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6, bl_type_7,  &
    ! OUT data required for tracer mixing :
      kent, we_lim, t_frac, zrzi, kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,  &
    ! OUT data required elsewhere in UM system :
      zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0,leonard_kl_tke                      &
      )

    if ( smagorinsky ) then

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
          weight1 = z_theta(i,1,k) - z_theta(i,1,k-1)
          weight2 = z_theta(i,1,k) - z_rho(i,1,k)
          weight3 = z_rho(i,1,k)   - z_theta(i,1,k-1)
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
        weight1 = z_theta(i,1,k) - z_theta(i,1,k-1)
        weight2 = z_theta(i,1,k) - z_rho(i,1,k)
        ! Assume rhokm(BL_LEVELS+1) is zero
        rhokm_mix(i,1,k) = (weight2/weight1) * rhokm(i,1,k)
        ! Scale exchange coefficients by 1/dz factor, as is done for
        ! rhokh_mix in bdy_impl4
        rhokm_mix(i,1,k) = rhokm_mix(i,1,k) * rdz_charney_grid(i,1,k)
      end do

      ! If full stress term is used, diagonal terms are doubled.
      if (fullstress) then
        do k = 1, bl_levels
          do i = 1, seg_len
            rhokm_mix(i,1,k) = 2.0 * rhokm_mix(i,1,k)
          end do
        end do
      end if

      kent_dummy = 2
      zeroes_2d = 0.0_r_um
      zeroes_ent = 0.0_r_um
      do k = 1, bl_levels
        do i = 1, seg_len
          w_mixed(i,1,k) = w(i,1,k)
        end do
      end do

      call  tr_mix (                                                           &
           ! IN fields
           r_theta_levels, r_rho_levels, pdims,                                &
           bl_levels, alpha_cd,                                                &
           rhokm_mix(1:seg_len,1:1,2:bl_levels),                               &
           rhokm_mix(1:seg_len,1:1,1),                                         &
           dtrdz_charney_grid, zeroes_2d, zeroes_2d, kent_dummy,               &
           zeroes_ent, zeroes_ent, zeroes_ent, kent_dummy,                     &
           zeroes_ent, zeroes_ent, zeroes_ent,                                 &
           zeroes_2d, zeroes_2d, real(z_rho,r_um),                             &
           ! INOUT / OUT fields
           w_mixed, w_flux, surf_dep_flux                                      &
           )

      do k = 1, bl_levels
        do i = 1, seg_len
          dw_bl(map_wth(1,i)+k) = w_mixed(i,1,k) - w(i,1,k)
        end do
      end do

    end if

    !----------------------------------------------------------------!
    ! Update variables
    !----------------------------------------------------------------!
    ! Set dummy zhnl value to avoid undefined reference in interpolation later
    do i = 1, seg_len
      surf_interp(map_surf(1,i)+5) = 0.0_r_def
    end do

    do k=1,bl_levels
      do i = 1, seg_len
        bq_bl(map_wth(1,i) + k-1) = bq_gb(i,1,k)
        bt_bl(map_wth(1,i) + k-1) = bt_gb(i,1,k)
        dtrdz_tq_bl(map_wth(1,i) + k) = dtrdz_charney_grid(i,1,k)
      end do
    end do

    ! Update variables for formdrag scheme
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

    ! Update variables for convection scheme
    do i = 1, seg_len
      zh_2d(map_2d(1,i))            = zh(i,1)
      ntml_2d(map_2d(1,i))          = ntml(i,1)
      if (cumulus(i,1)) then
        cumulus_2d(map_2d(1,i))     = 1_i_def
      else
        cumulus_2d(map_2d(1,i))     = 0_i_def
      end if
      z_lcl(map_2d(1,i))            = real(zlcl(i,1), r_def)
      inv_depth(map_2d(1,i))        = real(dzh(i,1), r_def)
      qcl_at_inv_top(map_2d(1,i))   = real(qcl_inv_top(i,1), r_def)
      if ( l_shallow(i,1) ) then
        shallow_flag(map_2d(1,i))   = 1_i_def
      else
        shallow_flag(map_2d(1,i))   = 0_i_def
      end if
      uw0_flux(map_2d(1,i))         = uw0(i,1)
      vw0_flux(map_2d(1,i))         = vw0(i,1)
      lcl_height(map_2d(1,i))       = zlcl_uv(i,1)
      parcel_top(map_2d(1,i))       = zhpar(i,1)
      level_parcel_top(map_2d(1,i)) = ntpar(i,1)
      wstar_2d(map_2d(1,i))         = wstar(i,1)
      thv_flux(map_2d(1,i))         = wthvs(i,1)
      parcel_buoyancy(map_2d(1,i))  = delthvu(i,1)
      qsat_at_lcl(map_2d(1,i))      = qsat_lcl(i,1)
    end do

    ! Update variables for boundary layer scheme
    do k = 1, bl_levels
      do i = 1, seg_len
        tke_bl(map_wth(1,i) + k-1) = tke_trb(i,1,k)
      end do
    end do

    if (bdy_tke /= bdy_tke_deardorff) then
      do k = 1, bl_levels
        do i = 1, seg_len
          tsq_bl(map_wth(1,i) + k-1) = tsq_trb(i,1,k)
          qsq_bl(map_wth(1,i) + k-1) = qsq_trb(i,1,k)
          cov_bl(map_wth(1,i) + k-1) = cov_trb(i,1,k)
        end do
      end do
      if (leonard_tke) then
        do k = 1, bl_levels
          do i = 1, seg_len
            leonard_klm_tke(map_wth(1,i) + k) = leonard_kl_tke(i,1,k,2)
            leonard_klh_tke(map_w3(1,i) + k-1) = leonard_kl_tke(i,1,k,1)
          end do
        end do
        do i = 1, seg_len
          leonard_klm_tke(map_wth(1,i)) = 0.0_r_def
          do k = bl_levels+1, nlayers
            leonard_klm_tke(map_wth(1,i) + k) = 0.0_r_def
            leonard_klh_tke(map_w3(1,i) + k-1) = 0.0_r_def
          end do
        end do
      end if
    end if

    if (shcu_buoy) then
      do i = 1, seg_len
        zhpar_shcu_2d(map_2d(1,i)) = zhpar_shcu(i,1)
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
      end do
    end do

    ! Counter gradient flux of momentum is defined at Wtheta level
    ! but tentatively store in W3 data to convert to W2 data with set_wind
    do k = 2, bl_levels
      do i = 1, seg_len
        ! true defined position is map_wth(1,i) +k-1
        rhogamu_w3(map_w3(1,i) + k-2) = rhogamu(i,1,k)
        rhogamv_w3(map_w3(1,i) + k-2) = rhogamv(i,1,k)
      end do
    end do
    do i = 1, seg_len
      do k = bl_levels, nlayers
        rhogamu_w3(map_w3(1,i) + k-1) = 0.0_r_def
        rhogamv_w3(map_w3(1,i) + k-1) = 0.0_r_def
      end do
    end do

    do i = 1, seg_len
      bl_type_ind(map_bl(1,i)+0) = bl_type_1(i,1)
      bl_type_ind(map_bl(1,i)+1) = bl_type_2(i,1)
      bl_type_ind(map_bl(1,i)+2) = bl_type_3(i,1)
      bl_type_ind(map_bl(1,i)+3) = bl_type_4(i,1)
      bl_type_ind(map_bl(1,i)+4) = bl_type_5(i,1)
      bl_type_ind(map_bl(1,i)+5) = bl_type_6(i,1)
      bl_type_ind(map_bl(1,i)+6) = bl_type_7(i,1)
    end do

    ! Update blended Smagorinsky diffusion coefficients
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
      do i = 1, seg_len
        do k = bl_levels+1, nlayers
          visc_m_blend(map_wth(1,i) + k) = 0.0_r_def
          visc_h_blend(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    endif

    if (.not. associated(rhogamu_bl, empty_real_data)) then
      do k = 2, bl_levels
        do i = 1, seg_len
          rhogamu_bl(map_wth(1,i) + k-1) = rhogamu(i,1,k)
        end do
      end do
      do i = 1, seg_len
        rhogamu_bl(map_wth(1,i)) = 0.0_r_def
        do k = bl_levels, nlayers
          rhogamu_bl(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(rhogamv_bl, empty_real_data)) then
      do k = 2, bl_levels
        do i = 1, seg_len
          rhogamv_bl(map_wth(1,i) + k-1) = rhogamv(i,1,k)
        end do
      end do
      do i = 1, seg_len
        rhogamv_bl(map_wth(1,i)) = 0.0_r_def
        do k = bl_levels, nlayers
          rhogamv_bl(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(rhogamt_bl, empty_real_data)) then
      do k = 2, bl_levels
        do i = 1, seg_len
          rhogamt_bl(map_w3(1,i) + k-1) = BL_diag%rhogamt(i,1,k)
        end do
      end do
      do i = 1, seg_len
        rhogamt_bl(map_w3(1,i)) = 0.0_r_def
        do k = bl_levels, nlayers
          rhogamt_bl(map_w3(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(rhogamq_bl, empty_real_data)) then
      do k = 2, bl_levels
        do i = 1, seg_len
          rhogamq_bl(map_w3(1,i) + k-1) = BL_diag%rhogamq(i,1,k)
        end do
      end do
      do i = 1, seg_len
        rhogamq_bl(map_w3(1,i)) = 0.0_r_def
        do k = bl_levels, nlayers
          rhogamq_bl(map_w3(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(tke_shr_prod, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          tke_shr_prod(map_wth(1,i) + k-1) = BL_diag%tke_shr_prod(i,1,k)
        end do
      end do
      do i = 1, seg_len
        tke_shr_prod(map_wth(1,i)) = 0.0_r_def
        do k = tke_levels, nlayers
          tke_shr_prod(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(tke_boy_prod, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          tke_boy_prod(map_wth(1,i) + k-1) = BL_diag%tke_boy_prod(i,1,k)
        end do
      end do
      do i = 1, seg_len
        tke_boy_prod(map_wth(1,i)) = 0.0_r_def
        do k = tke_levels, nlayers
          tke_boy_prod(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(tke_dissp, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          tke_dissp(map_wth(1,i) + k-1) = BL_diag%tke_dissp(i,1,k)
        end do
      end do
      do i = 1, seg_len
        tke_dissp(map_wth(1,i)) = 0.0_r_def
        do k = tke_levels, nlayers
          tke_dissp(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(master_length, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          master_length(map_wth(1,i) + k-1) = BL_diag%elm(i,1,k)
        end do
      end do
      do i = 1, seg_len
        ! Copy above vaule as mym_length does.
        master_length(map_wth(1,i)) = master_length(map_wth(1,i)+1)

        do k = tke_levels, nlayers
          master_length(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(sm25, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          sm25(map_wth(1,i) + k-1) = BL_diag%sm(i,1,k)
        end do
      end do
      do i = 1, seg_len
        sm25(map_wth(1,i)) = 0.0_r_def
        do k = tke_levels, nlayers
          sm25(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(sh25, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          sh25(map_wth(1,i) + k-1) = BL_diag%sh(i,1,k)
        end do
      end do
      do i = 1, seg_len
        sh25(map_wth(1,i)) = 0.0_r_def
        do k = tke_levels, nlayers
          sh25(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(dbdz, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          dbdz(map_wth(1,i) + k-1) = BL_diag%dbdz(i,1,k)
        end do
      end do
      do i = 1, seg_len
        dbdz(map_wth(1,i)) = 0.0_r_def
        do k = tke_levels, nlayers
          dbdz(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    if (.not. associated(dvdzm, empty_real_data)) then
      do k = 2, tke_levels
        do i = 1, seg_len
          dvdzm(map_wth(1,i) + k-1) = BL_diag%dvdzm(i,1,k)
        end do
      end do
      do i = 1, seg_len
        dvdzm(map_wth(1,i)) = 0.0_r_def
        do k = tke_levels, nlayers
          dvdzm(map_wth(1,i) + k) = 0.0_r_def
        end do
      end do
    end if

    ! deallocate diagnostics deallocated in atmos_physics2
    call dealloc_bl_expl(bl_diag)
    deallocate(BL_diag%gradrich)
    deallocate(land_index)
    deallocate(sd_orog)
    deallocate(ho2r2_orog_gb)
    deallocate(sil_orog_land_gb)

  end subroutine bl_hoc_exp_code

end module bl_hoc_exp_kernel_mod
