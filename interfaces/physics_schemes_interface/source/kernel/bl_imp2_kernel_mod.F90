!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the implicit UM boundary layer scheme.
!>
module bl_imp2_kernel_mod

  use argument_mod,           only : arg_type,                  &
                                     GH_FIELD, GH_SCALAR,       &
                                     GH_INTEGER, GH_REAL,       &
                                     GH_READ, GH_WRITE, GH_INC, &
                                     GH_READWRITE, DOMAIN,      &
                                     ANY_DISCONTINUOUS_SPACE_1, &
                                     ANY_DISCONTINUOUS_SPACE_2
  use section_choice_config_mod, only : cloud, cloud_um
  use blayer_config_mod,         only : fric_heating
  use cloud_config_mod,          only : scheme, scheme_smith, scheme_pc2, &
                                        scheme_bimodal,                   &
                                        i_bm_ez_opt, i_bm_ez_opt_entpar
  use constants_mod,             only : i_def, i_um, r_def, r_um, r_bl
  use fs_continuity_mod,         only : W3, Wtheta
  use kernel_mod,                only : kernel_type
  use timestepping_config_mod,   only : outer_iterations
  use mixing_config_mod,         only : leonard_term
  use physics_config_mod,        only : lowest_level,          &
                                        lowest_level_constant, &
                                        lowest_level_gradient, &
                                        lowest_level_flux
  use water_constants_mod,       only : lc, lf

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_imp2_kernel_type
    private
    type(arg_type) :: meta_args(57) = (/                                         &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! outer
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! loop
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA),                   &! dtheta_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! diss_u
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! diss_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dt_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_ci
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_s
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_area
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_ice
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_liq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cca
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ccw
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! rh_crit_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dsldzm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mix_len_bm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wvar
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! gradrinr
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! tau_dec_bm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! tau_hom_bm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! tau_mph_bm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rhokh_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! bq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! bt_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! moist_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! heat_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rdz_tq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! thetal_inc_leonard_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mt_inc_leonard_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! z_lcl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! inv_depth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! qcl_at_inv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zhsc_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! bl_type_ind
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! dqw_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! dtl_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! dqw_nt_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! dtl_nt_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! qw_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! tl_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! ct_ctq_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! fqw_star_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3)                        &! ftl_star_w3
         /)

    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: bl_imp2_code
  end type

  public :: bl_imp2_code

contains

  !> @brief Interface to the implicit UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     outer                Outer loop counter
  !> @param[in]     wetrho_in_wth        Wet density field in wth space
  !> @param[in]     exner_in_w3          Exner pressure field in density space
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     theta_star           Potential temperature after advection
  !> @param[in]     height_w3            Height of density space above surface
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in]     ntml_2d              Number of turbulently mixed levels
  !> @param[in]     cumulus_2d           Cumulus flag (true/false)
  !> @param[in,out] dtheta_bl            BL theta increment
  !> @param[in]     diss_u               Zonal Molecular dissipation rate
  !> @param[in]     diss_v               Meridional Molecular dissipation rate
  !> @param[in]     dt_conv              Convection temperature increment
  !> @param[in,out] m_v                  Vapour mixing ration after advection
  !> @param[in,out] m_cl                 Cloud liq mixing ratio after advection
  !> @param[in]     m_ci                 Cloud ice mixing ratio after advection
  !> @param[in,out] m_s                  Snow mixing ratio after advection
  !> @param[in,out] cf_area              Area cloud fraction
  !> @param[in,out] cf_ice               Ice cloud fraction
  !> @param[in,out] cf_liq               Liquid cloud fraction
  !> @param[in,out] cf_bulk              Bulk cloud fraction
  !> @param[in,out] cca                  Convective cloud amount (fraction)
  !> @param[in,out] ccw                  Convective cloud water (kg/kg) (can be ice or liquid)
  !> @param[in]     rh_crit_wth          Critical relative humidity
  !> @param[in]     dsldzm               Liquid potential temperature gradient in wth
  !> @param[in]     mix_len_bm           Turb length-scale for bimodal in wth
  !> @param[in]     wvar                 Vertical velocity variance in wth
  !> @param[in]     gradrinr             Gradient Richardson number in wth
  !> @param[in]     tau_dec_bm           Decorrelation time scale in wth
  !> @param[in]     tau_hom_bm           Homogenisation time scale in wth
  !> @param[in]     tau_mph_bm           Phase-relaxation time scale in wth
  !> @param[in]     rhokh_bl             Heat eddy diffusivity on BL levels
  !> @param[in]     bq_bl                Buoyancy parameter for moisture
  !> @param[in]     bt_bl                Buoyancy parameter for heat
  !> @param[in,out] moist_flux_bl        Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl         Vertical heat flux on BL levels
  !> @param[in]     dtrdz_tq_bl          dt/(rho*r*r*dz) in wth
  !> @param[in]     rdz_tq_bl            1/dz in w3
  !> @param[in]     thetal_inc_leonard_wth  Leonard term increment for water potential temperature
  !> @param[in]     mt_inc_leonard_wth   Leonard term increment for total moisture
  !> @param[in,out] z_lcl                Height of the LCL
  !> @param[in]     inv_depth            Depth of BL top inversion layer
  !> @param[in]     qcl_at_inv_top       Cloud water at top of inversion
  !> @param[in]     zh_nonloc            Depth of non-local BL scheme
  !> @param[in]     zh_2d                Total BL depth
  !> @param[in]     zhsc_2d              Height of decoupled layer top
  !> @param[in]     bl_type_ind          Diagnosed BL types
  !> @param[in,out] dqw_wth              Total increment to total water
  !> @param[in,out] dtl_wth              Total increment to liquid temperature
  !> @param[in]     dqw_nt_wth           Non-turbulent increment to total water
  !> @param[in]     dtl_nt_wth           Non-turbulent increment to liquid temperature
  !> @param[in,out] qw_wth               Total water
  !> @param[in,out] tl_wth               Liquid temperature
  !> @param[in]     ct_ctq_wth           Coefficient in predictor-corrector
  !> @param[in,out] fqw_star_w3          Total moisture flux predictor
  !> @param[in,out] ftl_star_w3          Heat flux predictor
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile             Number of DOFs per cell for tiles
  !> @param[in]     undf_tile            Number of total DOFs for tiles
  !> @param[in]     map_tile             Dofmap for cell for surface tiles
  !> @param[in]     ndf_pft              Number of DOFs per cell for PFTs
  !> @param[in]     undf_pft             Number of total DOFs for PFTs
  !> @param[in]     map_pft              Dofmap for cell for PFTs
  !> @param[in]     ndf_sice             Number of DOFs per cell for sice levels
  !> @param[in]     undf_sice            Number of total DOFs for sice levels
  !> @param[in]     map_sice             Dofmap for cell for sice levels
  !> @param[in]     ndf_soil             Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil            Number of total DOFs for soil levels
  !> @param[in]     map_soil             Dofmap for cell for soil levels
  !> @param[in]     ndf_smtile           Number of DOFs per cell for soil levels and tiles
  !> @param[in]     undf_smtile          Number of total DOFs for soil levels and tiles
  !> @param[in]     map_smtile           Dofmap for cell for soil levels and tiles
  !> @param[in]     ndf_bl               Number of DOFs per cell for BL types
  !> @param[in]     undf_bl              Number of total DOFs for BL types
  !> @param[in]     map_bl               Dofmap for cell for BL types
  subroutine bl_imp2_code(nlayers, seg_len,                   &
                          outer,                              &
                          loop,                               &
                          wetrho_in_wth,                      &
                          exner_in_w3,                        &
                          exner_in_wth,                       &
                          theta_star,                         &
                          height_w3,                          &
                          height_wth,                         &
                          ntml_2d,                            &
                          cumulus_2d,                         &
                          dtheta_bl,                          &
                          diss_u,                             &
                          diss_v,                             &
                          dt_conv,                            &
                          m_v,                                &
                          m_cl,                               &
                          m_ci,                               &
                          m_s,                                &
                          cf_area,                            &
                          cf_ice,                             &
                          cf_liq,                             &
                          cf_bulk,                            &
                          cca,                                &
                          ccw,                                &
                          rh_crit_wth,                        &
                          dsldzm,                             &
                          mix_len_bm,                         &
                          wvar,                               &
                          gradrinr,                           &
                          tau_dec_bm,                         &
                          tau_hom_bm,                         &
                          tau_mph_bm,                         &
                          rhokh_bl,                           &
                          bq_bl,                              &
                          bt_bl,                              &
                          moist_flux_bl,                      &
                          heat_flux_bl,                       &
                          dtrdz_tq_bl,                        &
                          rdz_tq_bl,                          &
                          thetal_inc_leonard_wth,             &
                          mt_inc_leonard_wth,                 &
                          z_lcl,                              &
                          inv_depth,                          &
                          qcl_at_inv_top,                     &
                          zh_nonloc,                          &
                          zh_2d,                              &
                          zhsc_2d,                            &
                          bl_type_ind,                        &
                          dqw_wth, dtl_wth, dqw_nt_wth,       &
                          dtl_nt_wth, qw_wth, tl_wth,         &
                          ct_ctq_wth, fqw_star_w3,            &
                          ftl_star_w3,                        &
                          ndf_wth,                            &
                          undf_wth,                           &
                          map_wth,                            &
                          ndf_w3,                             &
                          undf_w3,                            &
                          map_w3,                             &
                          ndf_2d,                             &
                          undf_2d,                            &
                          map_2d,                             &
                          ndf_bl, undf_bl, map_bl)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use atm_step_local, only: rhc_row_length, rhc_rows
    use bl_option_mod, only: l_noice_in_turb, puns, pstb, on
    use cloud_inputs_mod, only: i_cld_vn, i_pc2_init_logic, forced_cu, &
                                i_cld_area
    use cv_run_mod, only: l_param_conv
    use free_tracers_inputs_mod, only: l_wtrac
    use gen_phys_inputs_mod, only: l_mr_physics
    use mphys_inputs_mod, only: l_mcr_qcf2, l_casim
    use nlsizes_namelist_mod, only: bl_levels
    use pc2_constants_mod, only: i_cld_smith, i_cld_pc2,            &
                                 pc2init_logic_smooth, acf_off, i_cld_bimodal
    use planet_constants_mod, only: p_zero, kappa, planet_radius, g => g_bl,   &
                                    cp => cp_bl, lcrcp
    use timestep_mod, only: timestep

    ! subroutines used
    use bl_lsp_mod, only: bl_lsp
    use bl_diags_mod, only: strnewbldiag
    use bm_ctl_mod, only: bm_ctl
    use bdy_impl4_mod, only: bdy_impl4
    use pc2_bl_forced_cu_mod, only: pc2_bl_forced_cu
    use pc2_bl_inhom_ice_mod, only: pc2_bl_inhom_ice
    use pc2_delta_hom_turb_mod, only: pc2_delta_hom_turb
    use ls_arcld_mod, only: ls_arcld

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

    integer(kind=i_def), intent(in) :: ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl,seg_len)

    real(kind=r_def), dimension(undf_w3),  intent(inout) :: moist_flux_bl,     &
                                                            heat_flux_bl,      &
                                                            fqw_star_w3,       &
                                                            ftl_star_w3,       &
                                                            rhokh_bl
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dtheta_bl,         &
                                                            m_v, m_cl, m_s,    &
                                                            cf_area, cf_ice,   &
                                                            cf_liq, cf_bulk,   &
                                                            dqw_wth, dtl_wth,  &
                                                            qw_wth, tl_wth,    &
                                                            cca, ccw
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: exner_in_w3,        &
                                                           height_w3,          &
                                                           rdz_tq_bl,          &
                                                           diss_u, diss_v
    real(kind=r_def), dimension(undf_wth), intent(in) :: wetrho_in_wth,        &
                                                         exner_in_wth,         &
                                                         theta_star,           &
                                                         height_wth,           &
                                                         dqw_nt_wth,           &
                                                         dtl_nt_wth,           &
                                                         ct_ctq_wth,           &
                                                         dt_conv,              &
                                                         rh_crit_wth,          &
                                                         bq_bl, bt_bl,         &
                                                         dtrdz_tq_bl,          &
                                                         dsldzm,               &
                                                         mix_len_bm,           &
                                                         wvar, m_ci,            &
                                                         gradrinr,             &
                                                         tau_dec_bm,           &
                                                         tau_hom_bm,           &
                                                         tau_mph_bm,           &
                                                         mt_inc_leonard_wth,   &
                                                         thetal_inc_leonard_wth

    integer(kind=i_def), dimension(undf_2d), intent(in) :: ntml_2d,            &
                                                           cumulus_2d

    real(kind=r_def), dimension(undf_2d), intent(in) :: zh_2d,                &
                                                        zh_nonloc, zhsc_2d

    real(kind=r_def), intent(inout) :: z_lcl(undf_2d)
    real(kind=r_def), intent(in) :: inv_depth(undf_2d)
    real(kind=r_def), intent(in) :: qcl_at_inv_top(undf_2d)
    integer(kind=i_def), dimension(undf_bl), intent(in) :: bl_type_ind

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: k, i

    ! local switches and scalars
    integer(i_um) :: error_code
    logical :: l_correct

    type(strnewbldiag) :: bl_diag

    ! profile fields from level 1 upwards
    real(r_um), dimension(seg_len,1,nlayers) ::                              &
         p_rho_levels, z_theta, rhcpt, t_latest, q_latest, qcl_latest,       &
         qcf_latest, qcf2_latest, area_cloud_fraction,                       &
         cf_latest, cfl_latest, cff_latest, t_earliest, q_earliest,          &
         qcl_earliest, qcf_earliest, cf_earliest, cfl_earliest,              &
         cff_earliest, qt_force, tl_force, t_inc_pc2, q_inc_pc2, qcl_inc_pc2,&
         bcf_inc_pc2, cfl_inc_pc2, sskew, svar_turb, svar_bm, qcf_total,     &
         ri_bm, tgrad_in, mix_len_in, tau_dec_in, tau_hom_in, tau_mph_in,    &
         wvar_in
    real(r_bl), dimension(seg_len,1,nlayers) ::                              &
         r_rho_levels, rho_wet_tq

    ! Extra bimodal scheme diagnostics not yet implemented in LFRic
    real(r_um), dimension(1,1,1,1) :: entzone,                               &
                                      sl_modes, qw_modes, rh_modes, sd_modes
    ! Switch set to false to indicate not to calculate the above
    logical, parameter :: l_calc_diag = .FALSE.

    ! profile field on boundary layer levels
    real(r_bl), dimension(seg_len,1,bl_levels) :: fqw, ftl, rhokh,           &
         bq_gb, bt_gb, dtrdz_charney_grid, rdz_charney_grid, rhokh_mix, qw,  &
         tl, dqw, dtl, fqw_star, ftl_star

    ! profile fields on u/v points and all levels
    real(r_bl), dimension(seg_len,1,nlayers) :: r_u, r_v

    ! profile fields on u/v points and BL levels
    real(r_bl), dimension(seg_len,1,bl_levels) :: taux, tauy,                &
         rhokm_u, rhokm_v, dissip_u, dissip_v, taux_star, tauy_star, cq_cm_u,&
         cq_cm_v, ct_ctq, dqw_nt, dtl_nt, du_star, dv_star

    ! profile fields from level 2 upwards
    real(r_bl), dimension(seg_len,1,2:bl_levels) :: rdz_u, rdz_v

    ! profile fields from level 0 upwards
    real(r_um), dimension(seg_len,1,0:nlayers) ::                            &
         p_theta_levels, p_rho_minus_one
    real(r_bl), dimension(seg_len,1,0:nlayers) :: r_theta_levels

    ! single level real fields
    real(r_um), dimension(seg_len,1) ::                                      &
         zh, bl_type_3, bl_type_4, bl_type_6, bl_type_7, zhnl, zlcl_mix,     &
         zlcl, dzh, qcl_inv_top, zhsc
    real(r_bl), dimension(seg_len,1) :: gamma1, gamma2, fric_heating_blyr,   &
         fric_heating_incv

    ! single level integer fields
    integer(i_um), dimension(seg_len,1) :: ntml, nblyr

    ! single level logical fields
    logical, dimension(seg_len,1) :: cumulus

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    real(r_um), dimension(seg_len,1,nlayers) :: cca0, ccw0

    real(r_um), dimension(seg_len,1) :: xx_cos_theta_latitude

    integer(i_um), dimension(seg_len,1) :: lcbase0, ccb0, cct0

    ! parameters for new BL solver
    real(r_bl) :: pnonl,p1,p2
    real(r_bl), dimension(seg_len) :: i1, e1, e2
    real(r_bl), parameter :: sqrt2 = sqrt(2.0_r_bl)

    real(r_bl) :: weight1, weight2, weight3, ftl_m, fqw_m,     &
         f_buoy_m, dissip_mol, fric_heating_inc, z_blyr

    real(r_um), parameter :: qcl_max_factor = 0.1_r_um

    integer(i_um) :: large_levels
    integer(i_um), parameter :: levels_per_level = 3

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    ! Assuming map_wth(1) points to level 0
    ! and map_w3(1) points to level 1
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      do k = 1, nlayers
        ! height of rho levels from centre of planet
        r_rho_levels(i,1,k) = height_w3(map_w3(1,i) + k-1) + planet_radius
      end do
    end do
    if (loop == 1) then
      do i = 1, seg_len
        do k = 1, bl_levels
          dqw_nt(i,1,k) = dqw_nt_wth(map_wth(1,i) + k)
          dtl_nt(i,1,k) = dtl_nt_wth(map_wth(1,i) + k)
          dtrdz_charney_grid(i,1,k) = dtrdz_tq_bl(map_wth(1,i) + k)
        end do
      end do
    else
      do i = 1, seg_len
        do k = 2, bl_levels
          fqw_star(i,1,k) = fqw_star_w3(map_w3(1,i) + k-1)
          ftl_star(i,1,k) = ftl_star_w3(map_w3(1,i) + k-1)
        end do
      end do
    end if

    !-----------------------------------------------------------------------
    ! Things passed from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      do k = 1, bl_levels
        rhokh(i,1,k) = rhokh_bl(map_w3(1,i) + k-1)
        fqw(i,1,k) = moist_flux_bl(map_w3(1,i) + k-1)
        ftl(i,1,k) = heat_flux_bl(map_w3(1,i) + k-1)
        dqw(i,1,k) = dqw_wth(map_wth(1,i) + k)
        dtl(i,1,k) = dtl_wth(map_wth(1,i) + k)
        qw(i,1,k) = qw_wth(map_wth(1,i) + k)
        tl(i,1,k) = tl_wth(map_wth(1,i) + k)
        ct_ctq(i,1,k) = ct_ctq_wth(map_wth(1,i) + k)
        rdz_charney_grid(i,1,k) = rdz_tq_bl(map_w3(1,i) + k-1)
      end do
    end do

    !-----------------------------------------------------------------------
    ! Increments / fields with increments added
    !-----------------------------------------------------------------------
    if (loop == 2 .and. scheme == scheme_pc2) then
      do i = 1, seg_len
        do k = 1, nlayers
          t_earliest(i,1,k) = theta_star(map_wth(1,i) + k)   &
                            * exner_in_wth(map_wth(1,i) + k) &
                            + dt_conv(map_wth(1,i) + k)
          q_earliest(i,1,k) = m_v(map_wth(1,i) + k)
        end do
      end do
    end if
    !-----------------------------------------------------------------------
    ! Code from ni_imp_ctl
    !-----------------------------------------------------------------------
    bl_diag%l_ftl = .true.
    bl_diag%l_fqw = .true.

    do i = 1, seg_len
      p1=bl_type_ind(map_bl(1,i)+0)*pstb+(1.0_r_bl-bl_type_ind(map_bl(1,i)+0))*puns
      p2=bl_type_ind(map_bl(1,i)+1)*pstb+(1.0_r_bl-bl_type_ind(map_bl(1,i)+1))*puns
      pnonl=max(p1,p2)
      i1(i) = (1.0_r_bl+1.0_r_bl/sqrt2)*(1.0_r_bl+pnonl)
      e1(i) = (1.0_r_bl+1.0_r_bl/sqrt2)*( pnonl + (1.0_r_bl/sqrt2) + &
              sqrt(pnonl*(sqrt2-1.0_r_bl)+0.5_r_bl) )
      e2(i) = (1.0_r_bl+1.0_r_bl/sqrt2)*( pnonl+(1.0_r_bl/sqrt2) - &
              sqrt(pnonl*(sqrt2-1.0_r_bl)+0.5_r_bl))
      gamma1(i,1) = i1(i)
    end do

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

    call bdy_impl4 (                                                         &
         ! IN levels, switches
         bl_levels,  l_correct,                                              &
         ! IN data :
         gamma1, gamma2, rhokm_u, rhokm_v,                                   &
         rdz_charney_grid, r_rho_levels, dtrdz_charney_grid,rdz_u,rdz_v,     &
         ct_ctq,cq_cm_u,cq_cm_v,dqw_nt,dtl_nt,                               &
         ! INOUT data :
         qw,tl,fqw,ftl,taux,tauy,fqw_star,ftl_star,taux_star,tauy_star,      &
         r_u,r_v,du_star,dv_star,dqw,dtl,rhokh,bl_diag,                      &
         ! OUT data
         t_latest,q_latest,rhokh_mix                                         &
         )

    if (loop == 2) then

      do i = 1, seg_len
        do k = 0, nlayers
          ! height of theta levels from centre of planet
          r_theta_levels(i,1,k) = height_wth(map_wth(1,i) + k) + planet_radius
        end do
        zh(i,1) = zh_2d(map_2d(1,i))
      end do
      do k = 1, nlayers
        do i = 1, seg_len
          ! height of levels above surface
          z_theta(i,1,k) = r_theta_levels(i,1,k)-r_theta_levels(i,1,0)
        end do
      end do

      if (fric_heating) then

        do i = 1, seg_len
          do k = 1, bl_levels
            rho_wet_tq(i,1,k) = wetrho_in_wth(map_wth(1,i) + k)
            bq_gb(i,1,k) = bq_bl(map_wth(1,i) + k-1)
            bt_gb(i,1,k) = bt_bl(map_wth(1,i) + k-1)
            dissip_u(i,1,k) = diss_u(map_w3(1,i) + k-1)
            dissip_v(i,1,k) = diss_v(map_w3(1,i) + k-1)
          end do
        end do

        !-----------------------------------------------------------------------
        ! First, estimate molecular dissipation rate by assuming steady state
        ! subgrid KE, so that     dissip_mol = dissip_rke + f_buoy
        ! where f_buoy is the buoyancy flux interpolated to theta levels
        !-----------------------------------------------------------------------
        ! Then convert dissipation rate to heating rate,
        ! noting that dissipation rates are zero on BL_LEVELS,
        ! and redistribute within the boundary layer to account
        ! for lack of BL mixing of these increments
        !-----------------------------------------------------------------------
        do i = 1, seg_len
          nblyr(i,1) = 1
          k = 1

          weight1 = r_rho_levels(i,1,k+1) - r_theta_levels(i,1,0)
          weight2 = r_theta_levels(i,1,k) - r_theta_levels(i,1,0)
          weight3 = r_rho_levels(i,1,k+1) - r_theta_levels(i,1,k)
          ftl_m = weight2 * ftl(i,1,k+1) + weight3 * ftl(i,1,k)
          fqw_m = weight2 * fqw(i,1,k+1) + weight3 * fqw(i,1,k)
          f_buoy_m = g*( bt_gb(i,1,k)*(ftl_m/cp) +                             &
                         bq_gb(i,1,k)*fqw_m )/weight1

          dissip_mol = dissip_u(i,1,k)+dissip_v(i,1,k) + f_buoy_m
          fric_heating_inc = max (0.0_r_bl, timestep * dissip_mol             &
                                           / ( cp*rho_wet_tq(i,1,k) ) )

          ! Save level 1 heating increment for redistribution over
          ! boundary layer
          fric_heating_blyr(i,1) = fric_heating_inc *                          &
                       (r_rho_levels(i,1,2)-r_theta_levels(i,1,0))
        end do

        do k = 2, bl_levels-1
          do i = 1, seg_len
            weight1 = r_rho_levels(i,1,k+1) - r_rho_levels(i,1,k)
            weight2 = r_theta_levels(i,1,k) - r_rho_levels(i,1,k)
            weight3 = r_rho_levels(i,1,k+1) - r_theta_levels(i,1,k)
            ftl_m = weight2 * ftl(i,1,k+1) + weight3 * ftl(i,1,k)
            fqw_m = weight2 * fqw(i,1,k+1) + weight3 * fqw(i,1,k)

            f_buoy_m = g*( bt_gb(i,1,k)*(ftl_m/cp) +                           &
                           bq_gb(i,1,k)*fqw_m )/weight1

            dissip_mol = dissip_u(i,1,k)+dissip_v(i,1,k) + f_buoy_m
            fric_heating_incv(i,1) = max (0.0_r_bl, timestep * dissip_mol     &
                                                   / ( cp*rho_wet_tq(i,1,k) ) )

            if ( z_theta(i,1,k) <= zh(i,1) ) then
              !------------------------------------------------------
              ! Sum increments over boundary layer to avoid
              ! adding large increments in level 1
              !------------------------------------------------------
              nblyr(i,1) = k
              fric_heating_blyr(i,1) = fric_heating_blyr(i,1) +                &
                                       fric_heating_incv(i,1) *                &
                                  (r_rho_levels(i,1,k+1)-r_rho_levels(i,1,k))
            else
              t_latest(i,1,k) = t_latest(i,1,k) + fric_heating_incv(i,1)
              if (bl_diag%l_dtfric) then
                bl_diag%dtfric(i,1,k) = fric_heating_incv(i,1)
              end if
            end if

          end do
        end do

        !-----------------------------------------------------------------------
        ! Redistribute heating within the boundary layer to account
        ! for lack of BL mixing of these increments
        !-----------------------------------------------------------------------
        do i = 1, seg_len
          z_blyr = r_rho_levels(i,1,nblyr(i,1)+1)                              &
                 - r_theta_levels(i,1,0)
          fric_heating_blyr(i,1) = fric_heating_blyr(i,1) / z_blyr

          do k = 1, nblyr(i,1)

            ! Linearly decrease heating rate across surface layer
            fric_heating_inc = 2.0_r_bl * fric_heating_blyr(i,1) *            &
                              (1.0_r_bl-z_theta(i,1,k)/z_blyr)

            t_latest(i,1,k) = t_latest(i,1,k) + fric_heating_inc

            if (bl_diag%l_dtfric) then
              bl_diag%dtfric(i,1,k) = fric_heating_inc
            end if

          end do
        end do

      end if ! fric_heating

      !-----------------------------------------------------------------------
      ! Cloud section from ni_imp_ctl
      !-----------------------------------------------------------------------
      ! Create Tl and qT outside boundary layer levels
      do i = 1, seg_len
        do k = bl_levels+1, nlayers
          t_latest(i,1,k) = theta_star(map_wth(1,i) + k)   &
                            * exner_in_wth(map_wth(1,i) + k) &
                            + dt_conv(map_wth(1,i) + k)      &
                            - (lc * m_cl(map_wth(1,i) + k)) / cp
          q_latest(i,1,k) = m_v(map_wth(1,i) + k) + m_cl(map_wth(1,i) + k)
        end do
      end do

      if (cloud == cloud_um) then

        error_code=0

        do i = 1, seg_len
          cumulus(i,1) = (cumulus_2d(map_2d(1,i)) == 1_i_def)
          ntml(i,1) = ntml_2d(map_2d(1,i))
          dzh(i,1) = real(inv_depth(map_2d(1,i)), r_um)
          do k = 0, nlayers
            ! pressure on theta levels
            p_theta_levels(i,1,k) = p_zero*(exner_in_wth(map_wth(1,i) + k))**(1.0_r_def/kappa)
          end do
          do k = 1, nlayers
            qcl_latest(i,1,k) = m_cl(map_wth(1,i) + k)
            qcf_latest(i,1,k) = m_s(map_wth(1,i) + k)
            qcf2_latest(i,1,k) = m_ci(map_wth(1,i) + k)
          end do
        end do
        if (scheme == scheme_pc2) then
          do k = 1, nlayers
            do i = 1, seg_len
              qcl_earliest(i,1,k) = qcl_latest(i,1,k)
              qcf_earliest(i,1,k) = qcf_latest(i,1,k)
            end do
          end do
        end if

        do i = 1, seg_len
          do k = 1, nlayers
            ! Assign _latest with current updated values
            cfl_latest(i,1,k) = cf_liq(map_wth(1,i) + k)
            cff_latest(i,1,k) = cf_ice(map_wth(1,i) + k)
            cf_latest(i,1,k)  = cf_bulk(map_wth(1,i) + k)
          end do
        end do

        ! Remove qcf from T(liquid) and Q(vapour+liquid)
        if (.not. l_noice_in_turb)                                             &
             call bl_lsp( bl_levels, qcf_latest, q_latest, t_latest )

        ! Which cloud scheme are we using?
        ! 3 options are available:
        ! Specific switching off cloud scheme (i_cld_off)
        ! PC2 cloud scheme (i_cld_pc2)
        ! Otherwise use Smith cloud scheme (i_cld_smith)
        if ( i_cld_vn == i_cld_pc2) then

          do i = 1, seg_len
            zhnl(i,1) = zh_nonloc(map_2d(1,i))
            zlcl(i,1) = real(z_lcl(map_2d(1,i)), r_um)
            qcl_inv_top(i,1) = real(qcl_at_inv_top(map_2d(1,i)), r_um)
            bl_type_3(i,1) = bl_type_ind(map_bl(1,i)+2)
            bl_type_4(i,1) = bl_type_ind(map_bl(1,i)+3)
            bl_type_6(i,1) = bl_type_ind(map_bl(1,i)+5)
            zlcl_mix(i,1) = 0.0_r_um
          end do
          do k = 1, nlayers
            do i = 1, seg_len
              cf_earliest(i,1,k) = cf_latest(i,1,k)
              cfl_earliest(i,1,k) = cfl_latest(i,1,k)
              cff_earliest(i,1,k) = cff_latest(i,1,k)
            end do
          end do

          ! --------------------------------------------------------------------
          ! Inhomogeneous forcing of ice
          ! --------------------------------------------------------------------

          ! Calculate in-plume ice content (LS) by assuming that they are equal
          ! to the current values, except when the current value is not defined.
          ! Also calculate the forcing of ice content Q4F
          if (.not. l_casim) then
            !for casim there should be no mixing of ice and so latest
            !and earliest are the same
            call pc2_bl_inhom_ice( qcf_latest, qcf_earliest, cff_earliest,     &
                                   cff_latest, cf_latest)
          end if

          ! --------------------------------------------------------------------
          ! Homogeneous forcing of the liquid cloud
          ! --------------------------------------------------------------------

          ! Calculate forcing in qT and TL. Currently q_latest contains the
          ! vapour plus liquid content and q_earliest just the initial vapour
          ! content
          do k = 1, nlayers
            do i = 1, seg_len
              qt_force(i,1,k) = ( q_latest(i,1,k)                              &
                   - (q_earliest(i,1,k) + qcl_earliest(i,1,k)) )
              tl_force(i,1,k) = ( t_latest(i,1,k)                              &
                   - (t_earliest(i,1,k)- lc * qcl_earliest(i,1,k) / cp) )
            end do
          end do

          if ( leonard_term ) then
            ! Add turbulent increment due to the Leonard terms to the
            ! turbulent forcing of Tl and qw.  These increments were
            ! calculated earlier and get counted as
            ! part of the "non-turbulent" increment in imp_solver, but
            ! should contribute to the turbulent forcing of cloud here.
            do k = 1, bl_levels
              do i = 1, seg_len
                qt_force(i,1,k) = qt_force(i,1,k)                              &
                     + mt_inc_leonard_wth(map_wth(1,i) + k)
                tl_force(i,1,k) = tl_force(i,1,k)                              &
                     + thetal_inc_leonard_wth(map_wth(1,i) + k)                &
                     * exner_in_wth(map_wth(1,i) + k)
              end do
            end do
          end if

          ! Call homogeneous forcing routine
          call pc2_delta_hom_turb(                                             &
               ! INput variables
               p_theta_levels(1,1,1),                                          &
               ! INput variables
               t_earliest, q_earliest, qcl_earliest, cf_latest,                &
               cfl_latest, cff_latest, tl_force, qt_force,                     &
               ! OUTput variables
               t_inc_pc2, q_inc_pc2, qcl_inc_pc2, bcf_inc_pc2, cfl_inc_pc2,    &
               ! INput variables (other quantities)
               0.0_r_um, 0.0_r_um, l_mr_physics)

          if ( leonard_term ) then
            ! Leonard term increments were already added
            ! and are included in t_earliest / q_earliest here.
            ! We have added them to tl_force / qw_force above to include
            ! the cloud response to Leonard terms.
            ! But later the _latest variables are updated by adding _force
            ! to _earliest, so we'd be double-counting the Leonard terms.
            ! Therefore, subtract Leonard terms off the _force variables here:
            do k = 1, bl_levels
              do i = 1, seg_len
                qt_force(i,1,k) = qt_force(i,1,k)                              &
                     - mt_inc_leonard_wth(map_wth(1,i) + k)
                tl_force(i,1,k) = tl_force(i,1,k)                              &
                     - thetal_inc_leonard_wth(map_wth(1,i) + k)                &
                     * exner_in_wth(map_wth(1,i) + k)
              end do
            end do
          end if

          if ( i_pc2_init_logic/=pc2init_logic_smooth ) then
            ! Only do this removal of cloud at and below ntml if NOT using
            ! "smooth" PC2 initiation logic.  With "smooth" logic, we allow
            ! cloud to initiate below ntml, so also need to allow homogenous
            ! forcing to be applied to it for consistency.
            do k = 1, nlayers
              do i = 1, seg_len
              ! If below the LCL in cumulus or unstable but not shear-dominated
              ! boundary layers then we simply zero the liquid water content
              ! instead of using the homogeneous BL response. The forcings
              ! themselves are still applied but to q and T.
                if ( ( cumulus(i,1) .and. k  <=  ntml(i,1) ) .or.              &
                     ( forced_cu >= on .and. (bl_type_3(i,1) > 0.5_r_um        &
                     .or. bl_type_4(i,1) > 0.5_r_um )                          &
                     .and. z_theta(i,1,k)  <  zlcl(i,1) )  ) then
                  t_inc_pc2(i,1,k)   =  (-lcrcp) * qcl_earliest(i,1,k)
                  q_inc_pc2(i,1,k)   =  qcl_earliest(i,1,k)
                  qcl_inc_pc2(i,1,k) =  (-qcl_earliest(i,1,k))
                  cfl_inc_pc2(i,1,k) =  (-cfl_earliest(i,1,k))
                  bcf_inc_pc2(i,1,k) =  cff_latest(i,1,k)                      &
                       -cf_earliest(i,1,k)
                end if
              end do
            end do  ! k
          end if  ! ( i_pc2_init_logic/=pc2init_logic_smooth )

          ! To be consistent with the code above, set zlcl_mixed to
          ! prevent PC2 initiating cloud below this level
          ! Only needed for cumulus layer if l_param_conv is false - we can't
          ! initiate in cumulus when convection is active
          if (.not. l_param_conv) then
            do i = 1, seg_len
              if ( cumulus(i,1) ) zlcl_mix(i,1) = zlcl(i,1)
            end do
          end if

          ! With forced_cu, also do this for types 3 & 4
          if (forced_cu >= on) then
            do i = 1, seg_len
              if (bl_type_3(i,1) > 0.5_r_um .or. bl_type_4(i,1) > 0.5_r_um )   &
                   zlcl_mix(i,1) = zlcl(i,1)
            end do
          end if

          do k = 1, nlayers
            do i = 1, seg_len
              ! Update working version of temperature, moisture and cloud
              ! fields with increments from the PC2 homogeneous response.
              t_latest(i,1,k)   = t_earliest(i,1,k) + tl_force(i,1,k)          &
                   + t_inc_pc2(i,1,k)
              q_latest(i,1,k)   = q_earliest(i,1,k) + qt_force(i,1,k)          &
                   + q_inc_pc2(i,1,k)
              qcl_latest(i,1,k) = qcl_earliest(i,1,k)                          &
                   + qcl_inc_pc2(i,1,k)
              cfl_latest(i,1,k) = cfl_earliest(i,1,k)                          &
                   + cfl_inc_pc2(i,1,k)
              cf_latest(i,1,k)  =  cf_earliest(i,1,k)                          &
                   + bcf_inc_pc2(i,1,k)
              ! qcf_latest and cff_latest are not updated
            end do
          end do  ! k

          !------------------------------------------------------------
          ! Parametrize "forced cumulus clouds" at top of well-mixed BL
          !------------------------------------------------------------
          if (forced_cu >= on) then
            do i = 1, seg_len
              do k = 1, nlayers
                cca0(i,1,k) = cca(map_wth(1,i) + k)
                ccw0(i,1,k) = ccw(map_wth(1,i) + k)
              end do
            end do
            call pc2_bl_forced_cu( zhnl, dzh, zlcl, bl_type_3, bl_type_6,      &
                                   z_theta, qcl_inv_top,                       &
                                   cca0, ccw0, ccb0, cct0, lcbase0,            &
                                   cfl_latest, cf_latest,                      &
                                   qcl_latest, q_latest, t_latest, l_wtrac)
            do k = 1, nlayers
              do i = 1, seg_len
                cca(map_wth(1,i) + k) = cca0(i,1,k)
                ccw(map_wth(1,i) + k) = ccw0(i,1,k)
              end do
            end do
          end if  ! test on forced_cu

          ! --------------------------------------------------------------------
          ! Copy updated cloud fractions to the in/out variables
          ! --------------------------------------------------------------------
          do k = 1, nlayers
            do i = 1, seg_len
              ! For the moment set area cloud fraction
              ! to the bulk cloud fraction
              ! Ensure it has a value between 0.0 and 1.0
              area_cloud_fraction(i,1,k) = max(min(cf_latest(i,1,k),1.0_r_um),0.0_r_um)
            end do
          end do

          ! update BL prognostics
          if (outer == outer_iterations) then
            do i = 1, seg_len
              z_lcl(map_2d(1,i)) = zlcl_mix(i,1)
            end do
          end if !outer

        else if (i_cld_vn == i_cld_smith) then

          do i = 1, seg_len
            do k = 1, nlayers
              ! 3D RH_crit field
              rhcpt(i,1,k) = rh_crit_wth(map_wth(1,i) + k)
              p_rho_levels(i,1,k) = p_zero*(exner_in_w3(map_w3(1,i) + k-1))**(1.0_r_def/kappa)
            end do
          end do
          ! setup odd array which is on rho levels but without level 1
          do i = 1, seg_len
            p_rho_minus_one(i,1,0) = p_theta_levels(i,1,0)
          end do
          do k = 1, nlayers-1
            do i = 1, seg_len
              p_rho_minus_one(i,1,k) = p_rho_levels(i,1,k+1)
            end do
          end do
          do i = 1, seg_len
            p_rho_minus_one(i,1,nlayers) = 0.0_r_um
          end do

          ! --------------------------------------------------------------------
          ! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
          ! in boundary layer, calculate bulk_cloud fields from qT and qcf
          ! and calculate area_cloud fields.
          ! --------------------------------------------------------------------

          ! Determine number of sublevels for vertical gradient area cloud
          ! Want an odd number of sublevels per level: 3 is hardwired in do loop
          large_levels = ((nlayers - 2)*levels_per_level) + 2

          ! Add in qcf2 to diagnose cloud fraction for ice
          if (l_mcr_qcf2) then
            ! Two ice prognostics in use. Sum them together.
            do k = 1, nlayers
              do i = 1, seg_len
                qcf_total(i,1,k) = qcf_latest(i,1,k) + qcf2_latest(i,1,k)
              end do
            end do
          else ! l_mcr_qcf2
            ! Only one ice prognostic in use, so this is set to qcf_total
            do k = 1, nlayers
              do i = 1, seg_len
                qcf_total(i,1,k) = qcf_latest(i,1,k)
              end do
            end do
          end if ! l_mcr_qcf2

          call ls_arcld( p_theta_levels, rhcpt, p_rho_minus_one,               &
                         rhc_row_length, rhc_rows, bl_levels,                  &
                         levels_per_level, large_levels,                       &
                         xx_cos_theta_latitude,                                &
                         ntml, cumulus, l_mr_physics, qcf_total,               &
                         t_latest, q_latest, qcl_latest,                       &
                         area_cloud_fraction, cf_latest,                       &
                         cfl_latest, cff_latest,                               &
                         error_code)

        else if (i_cld_vn == i_cld_bimodal) then

          do i = 1, seg_len
            zhsc(i,1) = zhsc_2d(map_2d(1,i))
            bl_type_7(i,1) = bl_type_ind(map_bl(1,i)+6)
            do k = 1, nlayers
              tau_dec_in(i,1,k)  = tau_dec_bm(map_wth(1,i) + k)
              tau_hom_in(i,1,k)  = tau_hom_bm(map_wth(1,i) + k)
              tau_mph_in(i,1,k)  = tau_mph_bm(map_wth(1,i) + k)
              wvar_in(i,1,k)     = wvar(map_wth(1,i) + k )
            end do
          end do
          if (i_bm_ez_opt == i_bm_ez_opt_entpar) then
            ! Length-scale used for entraining parcel mode construction method
            do i = 1, seg_len
              do k = 1, nlayers
                mix_len_in(i,1,k)  = mix_len_bm(map_wth(1,i) + k)
              end do
            end do
          else
            ! SL-gradient used for stable-layer mode construction method
            do i = 1, seg_len
              do k = 1, nlayers
                tgrad_in(i,1,k)    = dsldzm(map_wth(1,i) + k)
              end do
            end do
          end if

          ! --------------------------------------------------------------------
          ! Section BL.4v Call bimodal cloud scheme to convert Tl and qT to T, q
          ! and qcl in boundary layer, calculate bulk_cloud fields from qT and
          ! qcf and calculate area_cloud fields.
          ! --------------------------------------------------------------------

          !---------------------------------------------------------------------
          ! Store Richardson number and components for zh_eff
          !---------------------------------------------------------------------
          do i = 1, seg_len
            do k = 1, bl_levels-1
              ! theta-level 1 is BL level 2, as BL has k=1 at the surface
              ri_bm(i,1,k) = gradrinr(map_wth(1,i) + k)
            end do
            do k = bl_levels, nlayers
              ! Set Ri-number to greater than 1 above bl_levels
              ri_bm(i,1,k) = 10.0_r_um
            end do
          end do

          if ( l_casim .and. l_mcr_qcf2 ) then
            ! Two ice prognostics in use. Sum them together
            do k = 1, nlayers
              do i = 1, seg_len
                qcf_total(i,1,k) = qcf_latest(i,1,k) + qcf2_latest(i,1,k)
              end do
            end do
          else ! l_casim / l_mcr_qcf2
            ! Only one ice prognostic in use, so this is set to qcf_total
            do k = 1, nlayers
              do i = 1, seg_len
                qcf_total(i,1,k) = qcf_latest(i,1,k)
              end do
            end do
          end if ! l_casim and l_mcr_qcf2

          call bm_ctl( p_theta_levels(:,:,1:), tgrad_in, wvar_in,              &
                       tau_dec_in, tau_hom_in, tau_mph_in, z_theta,            &
                       ri_bm, mix_len_in, zh, zhsc, dzh, bl_type_7,            &
                       nlayers, l_mr_physics, t_latest, cf_latest,             &
                       q_latest, qcf_total, qcl_latest,                        &
                       cfl_latest,cff_latest,                                  &
                       sskew, svar_turb, svar_bm, entzone,                     &
                       sl_modes, qw_modes, rh_modes, sd_modes,                 &
                       l_calc_diag, error_code )

          do k = 1, nlayers
            do i = 1, seg_len
              area_cloud_fraction(i,1,k)=cf_latest(i,1,k)
            end do
          end do

        end if

        ! update cloud fractions only if using cloud scheme
        do k = 1, nlayers
          do i = 1, seg_len
            cf_bulk(map_wth(1,i) + k) = cf_latest(i,1,k)
            cf_ice(map_wth(1,i) + k)  = cff_latest(i,1,k)
            cf_liq(map_wth(1,i) + k)  = cfl_latest(i,1,k)
            cf_area(map_wth(1,i) + k) = area_cloud_fraction(i,1,k)
          end do
        end do
        do i = 1, seg_len
          cf_ice(map_wth(1,i) + 0)  = cf_ice(map_wth(1,i) + 1)
          cf_liq(map_wth(1,i) + 0)  = cf_liq(map_wth(1,i) + 1)
          cf_bulk(map_wth(1,i) + 0) = cf_bulk(map_wth(1,i) + 1)
          cf_area(map_wth(1,i) + 0) = cf_area(map_wth(1,i) + 1)
        end do

      end if ! cloud_um

      !-----------------------------------------------------------------------
      ! Update main model prognostics
      !-----------------------------------------------------------------------
      do i = 1, seg_len
        do k = 1, nlayers
          ! potential temperature increment on theta levels
          dtheta_bl(map_wth(1,i) + k) = t_latest(i,1,k)                    &
                                      /  exner_in_wth(map_wth(1,i) + k)    &
                                      - theta_star(map_wth(1,i)+k)
          ! water vapour on theta levels
          m_v(map_wth(1,i) + k)  = q_latest(i,1,k)
          ! cloud liquid and ice water on theta levels
        end do
      end do
      if (cloud==cloud_um) then
        do k = 1, nlayers
          do i = 1, seg_len
            m_cl(map_wth(1,i) + k) = qcl_latest(i,1,k)
          end do
        end do
        do i = 1, seg_len
          m_cl(map_wth(1,i)) = m_cl(map_wth(1,i) + 1)
        end do
        if (.not. l_casim) then
          do k = 1, nlayers
            do i = 1, seg_len
              m_s(map_wth(1,i) + k) = qcf_latest(i,1,k)
            end do
          end do
          do i = 1, seg_len
            m_s(map_wth(1,i)) = m_s(map_wth(1,i) + 1)
          end do
        end if
      end if

      ! Update lowest-level values
      select case(lowest_level)
      case(lowest_level_constant)
        do i = 1, seg_len
          dtheta_bl(map_wth(1,i)) =                                            &
               t_latest(i,1,1) / exner_in_wth(map_wth(1,i) + 1)                &
               - theta_star(map_wth(1,i))
          m_v(map_wth(1,i))  = m_v(map_wth(1,i) + 1)
        end do

      case(lowest_level_gradient)
        do i = 1, seg_len
          dtheta_bl(map_wth(1,i)) =                                            &
               t_latest(i,1,1) / exner_in_wth(map_wth(1,i) + 1)                &
               - z_theta(i,1,1) * (                                            &
                   t_latest(i,1,2) / exner_in_wth(map_wth(1,i)+2)              &
                 - t_latest(i,1,1) / exner_in_wth(map_wth(1,i)+1)              &
                                  )                                            &
               / (z_theta(i,1,2) - z_theta(i,1,1)) - theta_star(map_wth(1,i))
          m_v(map_wth(1,i))  =                                                 &
               m_v(map_wth(1,i) + 1) - z_theta(i,1,1) * ( m_v(map_wth(1,i) + 2)&
             - m_v(map_wth(1,i) + 1) ) / (z_theta(i,1,2) - z_theta(i,1,1))
        end do

      case(lowest_level_flux)
        do i = 1, seg_len
          dtheta_bl(map_wth(1,i)) =                                            &
               t_latest(i,1,1) / exner_in_wth(map_wth(1,i) + 1)                &
               + ftl(i,1,1) / (cp * rhokh(i,1,1)) - theta_star(map_wth(1,i))
          m_v(map_wth(1,i))  = m_v(map_wth(1,i) + 1)                           &
               + fqw(i,1,1) / rhokh(i,1,1)
        end do
      end select

    else !loop=1

      do k = 1, bl_levels-1
        do i = 1, seg_len
          fqw_star_w3(map_w3(1,i) + k) = fqw_star(i,1,k+1)
          ftl_star_w3(map_w3(1,i) + k) = ftl_star(i,1,k+1)
        end do
      end do

    end if!loop

    do k = 0, bl_levels-1
      do i = 1, seg_len
        heat_flux_bl(map_w3(1,i)+k) = ftl(i,1,k+1)
        moist_flux_bl(map_w3(1,i)+k) = fqw(i,1,k+1)
        rhokh_bl(map_w3(1,i) + k) = rhokh(i,1,k+1)
      end do
    end do
    do k = 1, bl_levels
      do i = 1, seg_len
        dqw_wth(map_wth(1,i) + k) = dqw(i,1,k)
        dtl_wth(map_wth(1,i) + k) = dtl(i,1,k)
        qw_wth(map_wth(1,i) + k) = qw(i,1,k)
        tl_wth(map_wth(1,i) + k) = tl(i,1,k)
      end do
    end do

  end subroutine bl_imp2_code

end module bl_imp2_kernel_mod
