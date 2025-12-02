!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to CASIM microphysics scheme.

module casim_kernel_mod

use argument_mod,      only: arg_type,                  &
                             GH_FIELD, GH_REAL,         &
                             GH_READ, GH_WRITE,         &
                             GH_READWRITE,              &
                             ANY_DISCONTINUOUS_SPACE_1, &
                             ANY_DISCONTINUOUS_SPACE_2, &
                             CELL_COLUMN
use fs_continuity_mod, only: WTHETA, W3
use kernel_mod,        only: kernel_type
use empty_data_mod,    only: empty_real_data
use aerosol_config_mod, only: murk_prognostic
use microphysics_config_mod, only: casim_cdnc_opt, casim_cdnc_opt_external, &
                                   casim_cdnc_opt_fixed

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: casim_kernel_type
  private
  type(arg_type) :: meta_args(40) = (/                                      &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! ml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mi_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mr_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mg_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! ms_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! nl_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! nr_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! ni_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! ns_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! ng_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! w_phys
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! wetrho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! dry_rho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! height_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! height_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmv_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dml_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmi_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmr_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmg_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dms_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_rain_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_snow_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_graup_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! lsca_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! ls_rain_3d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! ls_snow_3d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! ls_graup_3d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! theta_inc
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! cloud_drop_no_conc
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! murk
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! refl_tot
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! refl_1km
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! superc_liq
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA)                        & ! superc_rain
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: casim_code
end type

public :: casim_code

contains

!> @brief Interface to the CASIM microphysics scheme
!>@details The CASIM Microphysics scheme calculates:
!>             1) Precipitation rates output to the surface
!>                and other physics schemes.
!>             2) Increments to the large scale prognostics
!>                due to cloud microphysical processes
!>                (e.g. latent heating and cooling).
!>         See UMDP50 for full scheme details
!> @param[in]     nlayers             Number of layers
!> @param[in]     mv_wth              Vapour mass mixing ratio
!> @param[in]     ml_wth              Liquid cloud mass mixing ratio
!> @param[in]     mi_wth              Ice cloud mass mixing ratio
!> @param[in]     mr_wth              Rain mass mixing ratio
!> @param[in]     mg_wth              Graupel mass mixing ratio
!> @param[in]     ms_wth              Snow mass mixing ratio
!> @param[in]     cfl_wth             Liquid cloud fraction
!> @param[in]     cff_wth             Ice cloud fraction
!> @param[in,out] nl_mphys            CASIM cloud-droplet number concentration
!> @param[in,out] nr_mphys            CASIM rain-drop number concentration
!> @param[in,out] ni_mphys            CASIM cloud-ice number concentration
!> @param[in,out] ns_mphys            CASIM snow number concentration
!> @param[in,out] ng_mphys            CASIM graupel number concentration
!> @param[in]     w_phys              'Vertical' wind in theta space
!> @param[in]     theta_in_wth        Potential temperature field
!> @param[in]     exner_in_wth        Exner pressure in potential temperature space
!> @param[in]     wetrho_in_w3        Wet density in density space
!> @param[in]     dry_rho_in_w3       Dry density in density space
!> @param[in]     height_w3           Height of density space levels above surface
!> @param[in]     height_wth          Height of theta levels above surface
!> @param[in,out] dmv_wth             Increment to vapour mass mixing ratio
!> @param[in,out] dml_wth             Increment to liquid cloud mass mixing ratio
!> @param[in,out] dmi_wth             Increment to ice cloud mass mixing ratio
!> @param[in,out] dmr_wth             Increment to rain mass mixing ratio
!> @param[in,out] dmg_wth             Increment to graupel mass mixing ratio
!> @param[in,out] dms_wth             Increment to snow mass mixing ratio
!> @param[in,out] ls_rain_2d          Large scale rain from twod_fields
!> @param[in,out] ls_snow_2d          Large scale snow from twod_fields
!> @param[in,out] ls_graup_2d         Large scale graupel from twod_fields
!> @param[in,out] lsca_2d             Large scale cloud amount (2d)
!> @param[in,out] ls_rain_3d          Large scale rain on model layers
!> @param[in,out] ls_snow_3d          Large scale snow on model layers
!> @param[in,out] ls_graup_3d         Large scale graupel on model layers
!> @param[in,out] theta_inc           Increment to theta
!> @param[in,out] cloud_drop_no_conc  In-cloud drop number for radiation
!> @param[in,out] refl_tot            Total radar reflectivity for diagnostic
!!                                     on all levels (dBZ)
!> @param[in,out] refl_1km            Radar reflectivity (dBZ) at 1km above the
!!                                     surface
!> @param[in,out] superc_liq          Supercooled liquid cloud mass mixing ratio
!> @param[in,out] superc_rain         Supercooled rain mass mixing ratio
!> @param[in]     ndf_wth             Number of degrees of freedom per cell for
!!                                     potential temperature space
!> @param[in]     undf_wth            Number unique of degrees of freedom for
!!                                     potential temperature space
!> @param[in]     map_wth             Dofmap for the cell at the base of the
!!                                     column for potential temperature space
!> @param[in]     ndf_w3              Number of degrees of freedom per cell for
!!                                     density space
!> @param[in]     undf_w3             Number unique of degrees of freedom for
!!                                     density space
!> @param[in]     map_w3              Dofmap for the cell at the base of the
!!                                     column for density space
!> @param[in]     ndf_2d              Number of degrees of freedom per cell for
!!                                     2D fields
!> @param[in]     undf_2d             Number unique of degrees of freedom for
!!                                     2D fields
!> @param[in]     map_2d              Dofmap for the cell at the base of the
!!                                     column for 2D fields

subroutine casim_code( nlayers,                     &
                       mv_wth,   ml_wth,   mi_wth,  &
                       mr_wth,   mg_wth,   ms_wth,  &
                       cfl_wth,  cff_wth,           &
                       nl_mphys, nr_mphys,          &
                       ni_mphys, ns_mphys, ng_mphys,&
                       w_phys,                      &
                       theta_in_wth,                &
                       exner_in_wth, wetrho_in_w3,  &
                       dry_rho_in_w3,               &
                       height_w3, height_wth,       &
                       dmv_wth,  dml_wth,  dmi_wth, &
                       dmr_wth,  dmg_wth,  dms_wth, &
                       ls_rain_2d, ls_snow_2d,      &
                       ls_graup_2d, lsca_2d,        &
                       ls_rain_3d, ls_snow_3d,      &
                       ls_graup_3d,                 &
                       theta_inc,                   &
                       cloud_drop_no_conc, murk,    &
                       refl_tot, refl_1km,          &
                       superc_liq, superc_rain,     &
                       ndf_wth, undf_wth, map_wth,  &
                       ndf_w3,  undf_w3,  map_w3,   &
                       ndf_2d,  undf_2d,  map_2d    )

    use constants_mod,              only: r_def, i_def, r_um, i_um
    use casim_diagnostics_mod,      only: ls_graup_3d_flag

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use timestep_mod,               only: timestep

    use atm_fields_bounds_mod,      only: pdims

    use planet_constants_mod,       only: p_zero, kappa, planet_radius
    use water_constants_mod,        only: tm

    use micro_main,                 only: shipway_microphysics
    use casim_switches,             only: its, ite, jts, jte, kts, kte, &
                                          ils, ile, jls, jle
    use generic_diagnostic_variables,                                  &
                                    only: allocate_diagnostic_space,   &
                                          deallocate_diagnostic_space, &
                                          casdiags
    use number_droplet_mod,         only: min_cdnc_sea_ice
    use mphys_air_density_mod,      only: mphys_air_density
    use mphys_radar_mod,            only: ref_lim
    use variable_precision,         only: wp
    use thresholds,            only: ql_tidy

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth,  ndf_w3,  ndf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_2d

    real(kind=r_def), intent(in),  dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mr_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mg_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ms_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: w_phys
    real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: wetrho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: dry_rho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: height_w3

    real(kind=r_def), intent(inout), dimension(undf_wth) :: nl_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: nr_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ni_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ns_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ng_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dml_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmi_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmr_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmg_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dms_wth
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_rain_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_snow_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_graup_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: lsca_2d
    real(kind=r_def), intent(inout), dimension(undf_wth) :: theta_inc
    real(kind=r_def), intent(inout), dimension(undf_wth) :: cloud_drop_no_conc
    real(kind=r_def), intent(inout), dimension(undf_wth) :: murk

    real(kind=r_def), intent(inout), dimension(undf_wth) :: ls_rain_3d
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ls_snow_3d

    real(kind=r_def), pointer, intent(inout) :: refl_tot(:)
    real(kind=r_def), pointer, intent(inout) :: refl_1km(:)

    real(kind=r_def), pointer, intent(inout) :: superc_liq(:)
    real(kind=r_def), pointer, intent(inout) :: superc_rain(:)
    real(kind=r_def), pointer, intent(inout) :: ls_graup_3d(:)

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    ! Local variables for the kernel
    real(wp), dimension(nlayers,1,1) ::                                        &
         qv_casim, qc_casim, qr_casim, nc_casim, nr_casim,                     &
         m3r_casim, qi_casim, qs_casim, qg_casim, ni_casim,                    &
         ns_casim, ng_casim, m3s_casim, m3g_casim,                             &
         th_casim,                                                             &
         aitken_sol_mass, aitken_sol_number, accum_sol_mass,                   &
         accum_sol_number, coarse_sol_mass, coarse_sol_number,                 &
         act_sol_liq_casim, act_sol_rain_casim, coarse_dust_mass,              &
         coarse_dust_number, act_insol_ice_casim,                              &
         act_sol_ice_casim, act_insol_liq_casim, accum_dust_mass,              &
         accum_dust_number, act_sol_number_casim,                              &
         act_insol_number_casim, aitken_sol_bk, accum_sol_bk,                  &
         coarse_sol_bk, pii_casim, p_casim,                                    &
         rho_casim, w_casim, tke_casim,                                        &
         dz_casim,                                                             &
         cfliq_casim, cfice_casim, cfsnow_casim,                               &
         cfrain_casim, cfgr_casim,                                             &
         dqv_casim, dqc_casim,  dqr_casim, dnc_casim,                          &
         dnr_casim, dm3r_casim, dqi_casim, dqs_casim,                          &
         dqg_casim, dni_casim, dns_casim,  dng_casim,                          &
         dm3s_casim, dm3g_casim, dth_casim,                                    &
         daitken_sol_mass, daitken_sol_number,                                 &
         daccum_sol_mass, daccum_sol_number,                                   &
         dcoarse_sol_mass, dcoarse_sol_number,                                 &
         dact_sol_liq_casim,   dact_sol_rain_casim,                            &
         dcoarse_dust_mass,    dcoarse_dust_number,                            &
         dact_insol_ice_casim, dact_sol_ice_casim,                             &
         dact_insol_liq_casim, daccum_dust_mass,                               &
         daccum_dust_number,   dact_sol_number_casim,                          &
         dact_insol_number_casim


    ! Local variables for the kernel
    real(r_um), parameter :: alt_1km = 1000.0_r_um ! metres

    real(r_um) :: t_work ! Local working temperature
    real(r_um) :: rrain, rsnow
    ! Scavenging rates, given in units of hr/mm
    real(r_um), parameter :: krain=2.0e-5_r_um, ksnow=2.0e-5_r_um

    logical :: l_refl_tot, l_refl_1km

    real(r_um), dimension(1,1,nlayers) ::                     &
         q_work, qcl_work, qcf_work, qrain_work, qcf2_work, qgraup_work,       &
         deltaz, rhodz_dry, rhodz_moist, rho_r2, dry_rho, r_rho_levels
    real(r_um), dimension(1,1,0:nlayers) :: r_theta_levels

    integer(i_um) :: k

    logical :: supercooled_layer(nlayers)

    !-------------------------------------------------------------------------
    ! End of Declarations
    !-------------------------------------------------------------------------

    ! Configure optional diagnostics
    casdiags % l_graupfall_3d = ls_graup_3d_flag

    ! Set CDNC for radiation here as we need the start of timestep value
    if (casim_cdnc_opt == casim_cdnc_opt_fixed) then
      do k = 0, nlayers
        if (cfl_wth(map_wth(1) + k) > 0.001_r_def) then
          cloud_drop_no_conc(map_wth(1) + k) = max(nl_mphys(map_wth(1) + k) / &
                                                   cfl_wth(map_wth(1) + k), &
                                                   min_cdnc_sea_ice)
        else
          cloud_drop_no_conc(map_wth(1) + k) = min_cdnc_sea_ice
        end if
      end do
    else if (casim_cdnc_opt == casim_cdnc_opt_external) then
      ! If we are getting the drop number from Glomap-clim or UKCA, then
      ! set the Casim drop number from this here. Ideally this would be done
      ! in casim_activate_kernel, but Glomap-clim and UKCA are both called
      ! after that has happened, hence why it needs to happen here.
      do k = 1, nlayers
        if (ml_wth(map_wth(1) + k) > ql_tidy) then
          nl_mphys( map_wth(1) + k) = cloud_drop_no_conc(map_wth(1) + k) * cfl_wth(map_wth(1) + k)
        else
          nl_mphys( map_wth(1) + k) = 0.0_r_def
        end if
      end do
    end if

    !-----------------------------------------------------------------------
    ! Initialisation of non-prognostic variables and arrays
    !-----------------------------------------------------------------------
    r_theta_levels(1,1,0) = height_wth(map_wth(1))+planet_radius
    do k = 1, nlayers
      ! height of levels from centre of planet
      r_rho_levels(1,1,k)   = height_w3(map_w3(1) + k-1) + planet_radius
      r_theta_levels(1,1,k) = height_wth(map_wth(1) + k) + planet_radius

      rho_r2(1,1,k)  = wetrho_in_w3(map_w3(1) + k-1) *                     &
                            ( r_rho_levels(1,1,k)**2 )
      dry_rho(1,1,k) = dry_rho_in_w3(map_w3(1) + k-1)
      ! Compulsory moist prognostics
      q_work(1,1,k)    = mv_wth(map_wth(1) + k)
      qcl_work(1,1,k)  = ml_wth(map_wth(1) + k)
      qcf_work(1,1,k)  = ms_wth(map_wth(1) + k)
      qcf2_work(1,1,k) = mi_wth(map_wth(1) + k)
      qrain_work(1,1,k) = mr_wth(map_wth(1) + k)
      qgraup_work(1,1,k) = mg_wth(map_wth(1) + k)
    end do     ! k

    ! calculate air density rhodz
    call mphys_air_density( r_theta_levels, r_rho_levels,                      &
                            dry_rho, rho_r2, pdims,                            &
                        q_work, qcl_work, qcf_work, qcf2_work,                 &
                        qrain_work, qgraup_work,                               &
                        rhodz_dry, rhodz_moist, deltaz )

    do k = 1, nlayers
      qv_casim(k,1,1) = mv_wth(map_wth(1) + k)
      qc_casim(k,1,1) = ml_wth(map_wth(1) + k)
      qr_casim(k,1,1) = mr_wth(map_wth(1) + k)
      nc_casim(k,1,1) = nl_mphys(map_wth(1) + k)
      nr_casim(k,1,1) = nr_mphys(map_wth(1) + k)
      m3r_casim(k,1,1) = 0.0_wp
      qi_casim(k,1,1) = mi_wth(map_wth(1) + k)
      qs_casim(k,1,1) = ms_wth(map_wth(1) + k)
      qg_casim(k,1,1) = mg_wth(map_wth(1) + k)
      ni_casim(k,1,1) = ni_mphys(map_wth(1) + k)
      ns_casim(k,1,1) = ns_mphys(map_wth(1) + k)
      ng_casim(k,1,1) = ng_mphys(map_wth(1) + k)
      m3s_casim(k,1,1) = 0.0_wp
      m3g_casim(k,1,1) = 0.0_wp
      th_casim(k,1,1) = theta_in_wth(map_wth(1) + k)
      aitken_sol_mass(k,1,1) = 0.0_wp
      aitken_sol_number(k,1,1) = 0.0_wp
      accum_sol_mass(k,1,1) =  0.0_wp
      accum_sol_number(k,1,1) = 0.0_wp
      coarse_sol_mass(k,1,1) = 0.0_wp
      coarse_sol_number(k,1,1) =  0.0_wp
      act_sol_liq_casim(k,1,1) = 0.0_wp
      act_sol_rain_casim(k,1,1) = 0.0_wp
      coarse_dust_mass(k,1,1) =  0.0_wp
      coarse_dust_number(k,1,1) = 0.0_wp
      act_insol_ice_casim(k,1,1) = 0.0_wp
      act_sol_ice_casim(k,1,1) = 0.0_wp
      act_insol_liq_casim(k,1,1) = 0.0_wp
      accum_dust_mass(k,1,1) =  0.0_wp
      accum_dust_number(k,1,1) =  0.0_wp
      act_sol_number_casim(k,1,1) =   0.0_wp
      act_insol_number_casim(k,1,1) = 0.0_wp
      aitken_sol_bk(k,1,1) = 0.0_wp
      accum_sol_bk(k,1,1) =   0.0_wp
      coarse_sol_bk(k,1,1) = 0.0_wp
      pii_casim(k,1,1) = exner_in_wth(map_wth(1) + k)
      p_casim(k,1,1) = p_zero*(exner_in_wth(map_wth(1) + k))               &
                                          **(1.0_wp/kappa)
      dz_casim(k,1,1)  = deltaz(1,1,k)
      rho_casim(k,1,1) = rhodz_dry(1,1,k) / dz_casim(k,1,1)
      w_casim(k,1,1) = w_phys(map_wth(1) + k)
      tke_casim(k,1,1) = 0.1_wp
      cfliq_casim(k,1,1) = cfl_wth(map_wth(1) + k)
      cfsnow_casim(k,1,1) = cff_wth(map_wth(1) + k)
      cfice_casim(k,1,1) = cfsnow_casim(k,1,1)

      dqv_casim(k,1,1) = 0.0_wp
      dqc_casim(k,1,1) = 0.0_wp
      dqr_casim(k,1,1) = 0.0_wp
      dnc_casim(k,1,1)  = 0.0_wp
      dnr_casim(k,1,1)  = 0.0_wp
      dm3r_casim(k,1,1) = 0.0_wp
      dqi_casim(k,1,1) = 0.0_wp
      dqs_casim(k,1,1) = 0.0_wp
      dqg_casim(k,1,1)  = 0.0_wp
      dni_casim(k,1,1) = 0.0_wp
      dns_casim(k,1,1)  = 0.0_wp
      dng_casim(k,1,1)  = 0.0_wp
      dm3s_casim(k,1,1) = 0.0_wp
      dm3g_casim(k,1,1) = 0.0_wp
      dth_casim(k,1,1) = 0.0_wp
      daitken_sol_mass(k,1,1) = 0.0_wp
      daitken_sol_number(k,1,1) = 0.0_wp
      daccum_sol_mass(k,1,1) = 0.0_wp
      daccum_sol_number(k,1,1) = 0.0_wp
      dcoarse_sol_mass(k,1,1) = 0.0_wp
      dcoarse_sol_number(k,1,1) = 0.0_wp
      dact_sol_liq_casim(k,1,1) = 0.0_wp
      dact_sol_rain_casim(k,1,1) = 0.0_wp
      dcoarse_dust_mass(k,1,1) = 0.0_wp
      dcoarse_dust_number(k,1,1) = 0.0_wp
      dact_insol_ice_casim(k,1,1) = 0.0_wp
      dact_sol_ice_casim(k,1,1) = 0.0_wp
      dact_insol_liq_casim(k,1,1) = 0.0_wp
      daccum_dust_mass(k,1,1) = 0.0_wp
      daccum_dust_number(k,1,1) = 0.0_wp
      dact_sol_number_casim(k,1,1) = 0.0_wp
      dact_insol_number_casim(k,1,1) = 0.0_wp
    end do     ! k

    cfrain_casim(nlayers,:,:)=0.0_wp
    cfgr_casim(nlayers,:,:)=0.0_wp
    do k =  nlayers-1, 1, -1
      !make cfrain the max of cfl in column
      cfrain_casim(k,1,1)=max(cfrain_casim(k+1,1,1),cfliq_casim(k,1,1),cfsnow_casim(k,1,1))
      !make graupel fraction
      cfgr_casim(k,1,1)=cfrain_casim(k,1,1)
    end do

    ! Set up diagnostic flags for CASIM
    l_refl_tot = .not. associated(refl_tot, empty_real_data)
    l_refl_1km = .not. associated(refl_1km, empty_real_data)

    if (l_refl_tot .or. l_refl_1km) casdiags % l_radar = .true.
    if (murk_prognostic) casdiags % l_snowfall_3d = .true.

    call allocate_diagnostic_space(its, ite, jts, jte, kts, kte)

    ! --------------------------------------------------------------------------
    ! this is the call to the CASIM microphysics
    ! Returns microphysical process rates
    ! --------------------------------------------------------------------------
    CALL shipway_microphysics( its, ite, jts, jte, kts, kte,  timestep,       &
                            qv_casim, qc_casim, qr_casim, nc_casim, nr_casim, &
                            m3r_casim, qi_casim, qs_casim, qg_casim, ni_casim,&
                            ns_casim, ng_casim, m3s_casim, m3g_casim,         &
                            th_casim,                                         &
                            aitken_sol_mass, aitken_sol_number,               &
                            accum_sol_mass,                                   &
                            accum_sol_number, coarse_sol_mass,                &
                            coarse_sol_number,                                &
                            act_sol_liq_casim, act_sol_rain_casim,            &
                            coarse_dust_mass,                                 &
                            coarse_dust_number, act_insol_ice_casim,          &
                            act_sol_ice_casim, act_insol_liq_casim,           &
                            accum_dust_mass,                                  &
                            accum_dust_number, act_sol_number_casim,          &
                            act_insol_number_casim, aitken_sol_bk,            &
                            accum_sol_bk,                                     &
                            coarse_sol_bk, pii_casim, p_casim,                &
                            rho_casim, w_casim, tke_casim,                    &
                            dz_casim,                                         &
                            cfliq_casim, cfice_casim, cfsnow_casim,           &
                            cfrain_casim, cfgr_casim,                         &
    !!                input variables above  || in/out variables below
                            dqv_casim, dqc_casim,  dqr_casim, dnc_casim,      &
                            dnr_casim, dm3r_casim, dqi_casim, dqs_casim,      &
                            dqg_casim, dni_casim, dns_casim,  dng_casim,      &
                            dm3s_casim, dm3g_casim, dth_casim,                &
                            daitken_sol_mass, daitken_sol_number,             &
                            daccum_sol_mass, daccum_sol_number,               &
                            dcoarse_sol_mass, dcoarse_sol_number,             &
                            dact_sol_liq_casim, dact_sol_rain_casim,          &
                            dcoarse_dust_mass,    dcoarse_dust_number,        &
                            dact_insol_ice_casim, dact_sol_ice_casim,         &
                            dact_insol_liq_casim, daccum_dust_mass,           &
                            daccum_dust_number,   dact_sol_number_casim,      &
                            dact_insol_number_casim,                          &
                            ils, ile,  jls, jle )

    ! Update murk for scavenging washout
    if (murk_prognostic) then
      ! Calculate scavenging rate in units of s/mm, to multiply by
      ! precip rate (mm/s)
      rrain = krain * timestep * 3600.0_r_um
      rsnow = ksnow * timestep * 3600.0_r_um
      do k = 1, nlayers
        murk(map_wth(1)+k) = murk(map_wth(1)+k) / &
             (1.0_r_um + rrain * casdiags%rainfall_3d(1,1,k) + &
                         rsnow * casdiags%snowfall_3d(1,1,k) )
      end do
      murk(map_wth(1)) = murk(map_wth(1)+1)
    end if

    ! CASIM Update theta and compulsory prognostic variables
    do k = 1, nlayers
      theta_inc(map_wth(1) + k) = dth_casim(k,1,1)
      dmv_wth(map_wth(1) + k ) = dqv_casim(k,1,1)
      dml_wth(map_wth(1) + k ) = dqc_casim(k,1,1)
      dmi_wth(map_wth(1) + k ) = dqi_casim(k,1,1)
      dms_wth(map_wth(1) + k ) = dqs_casim(k,1,1)
      dmr_wth( map_wth(1) + k) = dqr_casim(k,1,1)
      dmg_wth( map_wth(1) + k) = dqg_casim(k,1,1)
      nl_mphys( map_wth(1) + k) = nc_casim(k,1,1) +dnc_casim(k,1,1)
      nr_mphys( map_wth(1) + k) = nr_casim(k,1,1) +dnr_casim(k,1,1)
      ni_mphys( map_wth(1) + k) = ni_casim(k,1,1) +dni_casim(k,1,1)
      ns_mphys( map_wth(1) + k) = ns_casim(k,1,1) +dns_casim(k,1,1)
      ng_mphys( map_wth(1) + k) = ng_casim(k,1,1) +dng_casim(k,1,1)
    end do ! k (nlayers)

    ! Increment level 0 the same as level 1
    !  (as done in the UM)
    theta_inc(map_wth(1) + 0) = theta_inc(map_wth(1) + 1)
    dmv_wth(map_wth(1) + 0 ) = dmv_wth(map_wth(1) + 1 )
    dml_wth(map_wth(1) + 0 ) = dml_wth(map_wth(1) + 1 )
    dmi_wth(map_wth(1) + 0 ) = dmi_wth(map_wth(1) + 1 )
    dms_wth(map_wth(1) + 0 ) = dms_wth(map_wth(1) + 1 )
    dmr_wth( map_wth(1) + 0) = dmr_wth( map_wth(1) + 1)
    dmg_wth( map_wth(1) + 0) = dmg_wth( map_wth(1) + 1)
    nl_mphys( map_wth(1) + 0) = nl_mphys( map_wth(1) + 1)
    nr_mphys( map_wth(1) + 0) = nr_mphys( map_wth(1) + 1)
    ni_mphys( map_wth(1) + 0) = ni_mphys( map_wth(1) + 1)
    ns_mphys( map_wth(1) + 0) = ns_mphys( map_wth(1) + 1)
    ng_mphys( map_wth(1) + 0) = ng_mphys( map_wth(1) + 1)

    ! Copy ls_rain, ls_snow and ls_graup
    ls_rain_2d(map_2d(1))  = casdiags % SurfaceRainR(1,1)
    ls_snow_2d(map_2d(1))  = casdiags % SurfaceSnowR(1,1)
    ls_graup_2d(map_2d(1)) = casdiags % SurfaceGraupR(1,1)

    ! Copy 3D precipitation rate quantities
    do k = 1, nlayers
      ls_rain_3d(map_wth(1) + k)  = casdiags % rainfall_3d(1,1,k)
      ls_snow_3d(map_wth(1) + k)  = casdiags % snowonly_3d(1,1,k)
      if (ls_graup_3d_flag) then
        ls_graup_3d(map_wth(1) + k) = casdiags % graupfall_3d(1,1,k)
      end if
    end do

    ! Copy lsca_2d - like mphys_kernel_mod, use rain fraction
    ! from lowest model level
    lsca_2d(map_2d(1)) = cfrain_casim(1,1,1)

    if (l_refl_1km) then
      do k = 1, nlayers
        ! Select the first altitude above 1km (following what the UM does).
        if (height_wth(map_wth(1) + k) >= alt_1km ) then
          refl_1km(map_2d(1)) = casdiags % dbz_tot(1,1,k)
          exit
        end if
      end do
    end if

    if (l_refl_tot) then
      refl_tot(map_wth(1)) = ref_lim ! Set 0 level to -35 dBZ.
      do k = 1, nlayers
        refl_tot(map_wth(1) + k) = casdiags % dbz_tot(1,1,k)
      end do
    end if

    if (.not. associated(superc_liq, empty_real_data) .or.                     &
        .not. associated(superc_rain, empty_real_data) ) then
      do k = 1, nlayers
        t_work = exner_in_wth(map_wth(1) + k) * theta_in_wth(map_wth(1) + k)
        if (t_work < tm) then
          supercooled_layer(k) = .true.
        else
          supercooled_layer(k) = .false.
        end if
      end do

      if (.not. associated(superc_liq, empty_real_data) ) then
        do k = 1, nlayers
          if (supercooled_layer(k)) then
            superc_liq( map_wth(1) + k ) = ml_wth(map_wth(1) + k)
          else
            superc_liq( map_wth(1) + k ) = 0.0_r_um
          end if
        end do ! nlayers
      end if ! not assoc. superc_liq

      if (.not. associated(superc_rain, empty_real_data) ) then
        do k = 1, nlayers
          if (supercooled_layer(k)) then
            superc_rain( map_wth(1) + k ) = mr_wth(map_wth(1) + k)
          else
            superc_rain( map_wth(1) + k ) = 0.0_r_um
          end if
        end do ! nlayers
      end if ! not assoc. superc_rain
    end if ! not assoc. either superc species

    ! CASIM deallocate diagnostics
    call deallocate_diagnostic_space()
    ! (The above subroutine call sets casdiags % l_radar = .false.)

end subroutine casim_code

end module casim_kernel_mod
