!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Run the CFMIP Observation Simulator Package (COSP)
module cosp_kernel_mod

use argument_mod,      only : arg_type, &
                              GH_FIELD, GH_SCALAR, &
                              GH_REAL, GH_INTEGER, &
                              GH_READ, GH_WRITE, &
                              DOMAIN, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5
use fs_continuity_mod, only : W3, Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

type, public, extends(kernel_type) :: cosp_kernel_type
  private
  type(arg_type) :: meta_args(50) = (/ &
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! pressure_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! temperature_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! rho_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! height_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        & ! pressure_in_w3
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        & ! height_w3
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! mv
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! mcl
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! mcf
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! n_ice
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! radiative_cloud_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! frozen_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! radiative_conv_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! conv_liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! conv_frozen_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! conv_liquid_mmr
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! conv_frozen_mmr
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! conv_frozen_number
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! sigma_ml
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! sigma_mi
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! cloud_drop_no_conc
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! ls_rain_3d
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! conv_rain_3d
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! conv_snow_3d
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction
    arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! rand_seed
    arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! n_cloud_layer
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             & ! n_subcol_gen
    arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             & ! x1r
    arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             & ! x2r
    arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             & ! x1g
    arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             & ! x2g
    arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             & ! x4g
    ! Diagnostics (section cosp__)
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sunlit_mask
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! isccp_ctp_tau
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! calipso_low_cloud
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! calipso_low_cloud_mask
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! calipso_mid_cloud
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! calipso_mid_cloud_mask
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! calipso_high_cloud
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! calipso_high_cloud_mask
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), & ! calipso_cf_40_lvls_liq
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), & ! calipso_cf_40_lvls_ice
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), & ! calipso_cf_40_lvls_undet
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), & ! calipso_cf_40_lvls_mask
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), & ! calipso_total_backscatter
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! calipso_cfad_sr_40_lvls
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, Wtheta),                    & ! cloud_thermal_absorptivity
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, Wtheta)                     & ! cloud_solar_extinction
    /)

  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: cosp_code

end type cosp_kernel_type

public :: cosp_code

contains

!> @param[in]    nlayers                    Number of layers in Wtheta field
!> @param[in]    n_profile                  Number of columns
!> @param[in]    pressure_in_wth            Pressure in Wtheta space
!> @param[in]    temperature_in_wth         Temperature in Wtheta space
!> @param[in]    rho_in_wth                 Density in Wtheta space
!> @param[in]    height_wth                 Height in Wtheta space
!> @param[in]    pressure_in_w3             Pressure in W3 space
!> @param[in]    height_w3                  Height in W3 space
!> @param[in]    mv                         Water vapour field
!> @param[in]    mcl                        Large scale cloud liquid gridbox MMR
!> @param[in]    mcf                        Large scale cloud frozen gridbox MMR
!> @param[in]    radiative_cloud_fraction   Large scale cloud fraction
!> @param[in]    liquid_fraction            Large scale liquid cloud fraction
!> @param[in]    frozen_fraction            Large scale frozen cloud fraction
!> @param[in]    radiative_conv_fraction    Convective cloud fraction
!> @param[in]    conv_liquid_fraction       Convective liquid cloud fraction
!> @param[in]    conv_frozen_fraction       Convective frozen cloud fraction
!> @param[in]    conv_liquid_mmr            Convective liquid gridbox MMR
!> @param[in]    conv_frozen_mmr            Convective frozen gridbox MMR
!> @param[in]    conv_frozen_number         Convective frozen number conc
!> @param[in]    sigma_ml                   Fractional standard deviation of liquid condensate
!> @param[in]    sigma_mi                   Fractional standard deviation of ice condensate
!> @param[in]    cloud_drop_no_conc         Cloud Droplet Number Concentration
!> @param[in]    ls_rain_3d                 Large scale rain
!> @param[in]    conv_rain_3d               Convective rain
!> @param[in]    conv_snow_3d               Convective snow
!> @param[in]    lit_fraction               Lit fraction of the timestep
!> @param[in]    rand_seed                  Random seed field for cloud generator
!> @param[in]    n_cloud_layer              Number of cloud layers
!> @param[in]    n_subcol_gen               Number of cloud subcolumns
!> @param[in]    x1r                        Rain parameter
!> @param[in]    x2r                        Rain parameter
!> @param[in]    x1g                        Graupel parameter
!> @param[in]    x2g                        Graupel parameter
!> @param[in]    x4g                        Graupel parameter
!> @param[inout] sunlit_mask                COSP diagnostic
!> @param[inout] isccp_ctp_tau              COSP diagnostic
!> @param[inout] calipso_low_cloud          COSP diagnostic
!> @param[inout] calipso_low_cloud_mask     COSP diagnostic
!> @param[inout] calipso_mid_cloud          COSP diagnostic
!> @param[inout] calipso_mid_cloud_mask     COSP diagnostic
!> @param[inout] calipso_high_cloud         COSP diagnostic
!> @param[inout] calipso_high_cloud_mask    COSP diagnostic
!> @param[inout] calipso_cf_40_lvls_liq     COSP diagnostic
!> @param[inout] calipso_cf_40_lvls_ice     COSP diagnostic
!> @param[inout] calipso_cf_40_lvls_undet   COSP diagnostic
!> @param[inout] calipso_cf_40_lvls_mask    COSP diagnostic
!> @param[inout] calipso_total_backscatter  COSP diagnostic
!> @param[inout] calipso_cfad_sr_40_lvls    COSP diagnostic
!> @param[inout] cloud_thermal_absorptivity COSP diagnostic
!> @param[inout] cloud_solar_extinction     COSP diagnostic
!> @param[in]    ndf_wtheta                 No. DOFs per cell for Wtheta space
!> @param[in]    undf_wtheta                No. unique DOFs for Wtheta space
!> @param[in]    map_wtheta                 Dofmap for Wtheta column base cell
!> @param[in]    ndf_w3                     No. DOFs per cell for W3 space
!> @param[in]    undf_w3                    No. unique DOFs for W3 space
!> @param[in]    map_w3                     Dofmap for W3 column base cell
!> @param[in]    ndf_2d                     No. DOFs per cell for 2d space
!> @param[in]    undf_2d                    No. unique DOFs for 2d space
!> @param[in]    map_2d                     Dofmap for 2d column base cell
!> @param[in]    ndf_ptau                   No. DOFs per cell for ptau space
!> @param[in]    undf_ptau                  No. unique DOFs for ptau space
!> @param[in]    map_ptau                   Dofmap for ptau column base cell
!> @param[in]    ndf_40                     No. DOFs per cell for 40 space
!> @param[in]    undf_40                    No. unique DOFs for 40 space
!> @param[in]    map_40                     Dofmap for 40 column base cell
!> @param[in]    ndf_subcol                 No. DOFs per cell for subcol space
!> @param[in]    undf_subcol                No. unique DOFs for subcol space
!> @param[in]    map_subcol                 Dofmap for subcol column base cell
!> @param[in]    ndf_atb40                  No. DOFs per cell for atb40 space
!> @param[in]    undf_atb40                 No. unique DOFs for atb40 space
!> @param[in]    map_atb40                  Dofmap for atb40 column base cell
subroutine cosp_code(nlayers, n_profile, &
                     pressure_in_wth, temperature_in_wth, rho_in_wth, height_wth, &
                     pressure_in_w3, height_w3, &
                     mv, mcl, mcf, n_ice, &
                     radiative_cloud_fraction, &
                     liquid_fraction, frozen_fraction, &
                     radiative_conv_fraction, &
                     conv_liquid_fraction, conv_frozen_fraction, &
                     conv_liquid_mmr, conv_frozen_mmr, conv_frozen_number, &
                     sigma_ml, sigma_mi, cloud_drop_no_conc, &
                     ls_rain_3d, conv_rain_3d, conv_snow_3d, &
                     lit_fraction, rand_seed, n_cloud_layer, &
                     n_subcol_gen, x1r, x2r, x1g, x2g, x4g, &
                     sunlit_mask, &
                     isccp_ctp_tau, &
                     calipso_low_cloud, calipso_low_cloud_mask, &
                     calipso_mid_cloud, calipso_mid_cloud_mask, &
                     calipso_high_cloud, calipso_high_cloud_mask, &
                     calipso_cf_40_lvls_liq, &
                     calipso_cf_40_lvls_ice, &
                     calipso_cf_40_lvls_undet, &
                     calipso_cf_40_lvls_mask, &
                     calipso_total_backscatter, &
                     calipso_cfad_sr_40_lvls, &
                     cloud_thermal_absorptivity, &
                     cloud_solar_extinction, &
                     ndf_wtheta, undf_wtheta, map_wtheta, &
                     ndf_w3, undf_w3, map_w3, &
                     ndf_2d, undf_2d, map_2d, &
                     ndf_ptau, undf_ptau, map_ptau, &
                     ndf_40, undf_40, map_40, &
                     ndf_subcol, undf_subcol, map_subcol, &
                     ndf_atb40, undf_atb40, map_atb40)

  use empty_data_mod, only: empty_real_data
  use radiation_config_mod, only: i_cloud_ice_type_lw, i_cloud_liq_type_lw, &
                                  i_cloud_ice_type_sw, i_cloud_liq_type_sw, &
                                  cloud_vertical_decorr, &
                                  constant_droplet_effective_radius, &
                                  liu_aparam, liu_bparam
  use socrates_init_mod, only: i_cloud_representation, i_overlap, i_drop_re
  use socrates_runes, only: runes, StrDiag, &
                            ip_source_thermal, ip_source_illuminate, &
                            ip_inhom_mcica, ip_inhom_homogeneous
  use cosp_mod, only: cosp, &
                      n_cloudsat_levels, n_backscatter_bins,  &
                      n_isccp_tau_bins, n_isccp_pressure_bins
  use cosp_def_diag, only: CospDiag

  implicit none

  ! Dummy arguments
  integer(i_def), intent(in) :: nlayers, n_profile, n_subcol_gen
  real(r_def), intent(in) :: x1r, x2r, x1g, x2g, x4g

  integer(i_def), intent(in) :: ndf_wtheta, undf_wtheta
  integer(i_def), intent(in), dimension(ndf_wtheta, n_profile) :: map_wtheta
  real(r_def), intent(in), dimension(undf_wtheta), target :: &
    pressure_in_wth, temperature_in_wth, rho_in_wth, height_wth, &
    mv, mcl, mcf, n_ice, &
    radiative_cloud_fraction, liquid_fraction, frozen_fraction, &
    radiative_conv_fraction, conv_liquid_fraction, conv_frozen_fraction, &
    conv_liquid_mmr, conv_frozen_mmr, sigma_ml, sigma_mi, cloud_drop_no_conc, &
    ls_rain_3d, conv_rain_3d, conv_snow_3d, conv_frozen_number

  integer(i_def), intent(in) :: ndf_w3, undf_w3
  integer(i_def), intent(in), dimension(ndf_w3, n_profile) :: map_w3
  real(r_def), intent(in), dimension(undf_w3), target :: &
    pressure_in_w3, height_w3

  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in), dimension(ndf_2d, n_profile) :: map_2d
  integer(i_def), intent(in), dimension(undf_2d) :: rand_seed, n_cloud_layer
  real(r_def), intent(in), dimension(undf_2d) :: lit_fraction

  integer(i_def), intent(in) :: ndf_ptau, undf_ptau
  integer(i_def), intent(in), dimension(ndf_ptau, n_profile) :: map_ptau

  integer(i_def), intent(in) :: ndf_40, undf_40
  integer(i_def), intent(in), dimension(ndf_40, n_profile) :: map_40

  integer(i_def), intent(in) :: ndf_subcol, undf_subcol
  integer(i_def), intent(in), dimension(ndf_subcol, n_profile) :: map_subcol

  integer(i_def), intent(in) :: ndf_atb40, undf_atb40
  integer(i_def), intent(in), dimension(ndf_atb40, n_profile) :: map_atb40

  ! Conditional Diagnostics
  real(r_def), pointer, dimension(:), intent(inout) :: & ! 2d
    sunlit_mask, &
    calipso_low_cloud, calipso_low_cloud_mask, &
    calipso_mid_cloud, calipso_mid_cloud_mask, &
    calipso_high_cloud, calipso_high_cloud_mask
  real(r_def), pointer, dimension(:), intent(inout) :: & ! ptau
    isccp_ctp_tau
  real(r_def), pointer, dimension(:), intent(inout) :: & ! 40
    calipso_cf_40_lvls_liq, calipso_cf_40_lvls_ice, &
    calipso_cf_40_lvls_undet, calipso_cf_40_lvls_mask
  real(r_def), pointer, dimension(:), intent(inout) :: & ! subcol
    calipso_total_backscatter
  real(r_def), pointer, dimension(:), intent(inout) :: & ! atb40
    calipso_cfad_sr_40_lvls
  real(r_def), pointer, dimension(:), intent(inout) :: & ! wtheta
    cloud_thermal_absorptivity, cloud_solar_extinction

  ! Local variables for the kernel
  integer(i_def) :: n_profile_list, n_subcol_cosp
  integer(i_def), allocatable :: profile_list(:)
  integer(i_def) :: k, l
  integer(i_def) :: wth_0, wth_1, wth_last
  integer(i_def) :: w3_1, w3_last
  integer(i_def) :: twod_1, twod_last
  integer(i_def) :: ptau_1, ptau_last
  integer(i_def) :: c40_1, c40_last
  integer(i_def) :: subcol_0, subcol_last
  integer(i_def) :: atb40_1, atb40_last
  type(StrDiag) :: lw_diag, sw_diag
  real(r_def), pointer, dimension(:, :) :: p_full_levels
  real(r_def), pointer, dimension(:, :) :: p_half_levels
  real(r_def), pointer, dimension(:, :) :: hgt_full_levels
  real(r_def), pointer, dimension(:, :) :: hgt_half_levels
  real(r_def), pointer, dimension(:, :) :: density_full_levels
  real(r_def), pointer, dimension(:, :) :: d_mass
  real(r_def), pointer, dimension(:, :) :: t_full_levels
  real(r_def), pointer, dimension(:, :) :: q_full_levels
  real(r_def), pointer, dimension(:, :) :: cosp_crain
  real(r_def), pointer, dimension(:, :) :: cosp_csnow
  real(r_def), pointer, dimension(:) :: t_surf, p_surf, hgt_surf
  real(r_def), pointer, dimension(:) :: cosp_sunlit
  type(CospDiag) :: cosp_diag


  ! Set indexing
  wth_0 = map_wtheta(1,1)
  wth_1 = map_wtheta(1,1)+1
  wth_last = map_wtheta(1,1)+n_profile*(nlayers+1)-1
  w3_1 = map_w3(1,1)
  w3_last = map_w3(1,1)+n_profile*nlayers-1
  twod_1 = map_2d(1,1)
  twod_last = map_2d(1,1)+n_profile-1
  ptau_1 = map_ptau(1,1)
  ptau_last = map_ptau(1,1)+n_profile*n_isccp_tau_bins*n_isccp_pressure_bins-1
  c40_1 = map_40(1,1)
  c40_last = map_40(1,1)+n_profile*n_cloudsat_levels-1
  subcol_0 = map_subcol(1,1)
  subcol_last = map_subcol(1,1)+n_profile*(nlayers+1)*n_subcol_gen-1
  atb40_1 = map_atb40(1,1)
  atb40_last = map_atb40(1,1)+n_profile*n_cloudsat_levels*n_backscatter_bins-1

  ! Set pointers for LW diagnostic output from runes call
  allocate( &
    lw_diag%n_subcol_cloud(1:n_profile), &
    lw_diag%liq_subcol_scaling(1:nlayers, 1:n_subcol_gen, 1:n_profile), &
    lw_diag%total_cloud_fraction(1:nlayers, 1:n_profile), &
    lw_diag%liq_dim(1:nlayers, 1:n_profile), &
    lw_diag%liq_part_frac(1:nlayers, 1:n_profile), &
    lw_diag%liq_incloud_mmr(1:nlayers, 1:n_profile), &
    lw_diag%ice_re(1:nlayers, 1:n_profile), &
    lw_diag%ice_part_frac(1:nlayers, 1:n_profile), &
    lw_diag%ice_incloud_mmr(1:nlayers, 1:n_profile) )

  ! Diagnose cloud absorptivity at 11 micron
  lw_diag%cloud_absorptivity_wavelength = 11.0e-6_r_def
  if (.not. associated(cloud_thermal_absorptivity, empty_real_data)) then
    lw_diag%cloud_thermal_absorptivity(0:nlayers, 1:n_profile) &
      => cloud_thermal_absorptivity(wth_0:wth_last)
  else
    allocate( lw_diag%cloud_thermal_absorptivity(0:nlayers, 1:n_profile) )
  end if

  ! Run the Socrates cloud generator and calculate cloud diagnostics from LW
  do k=0, nlayers
    profile_list = pack( [(l, l=1, n_profile)], &
                         n_cloud_layer(twod_1:twod_last) == k )
    n_profile_list = size(profile_list)
    if (n_profile_list > 0) then
      call runes(n_profile_list, nlayers, lw_diag,                           &
        spectrum_name          = 'lw',                                       &
        i_source               = ip_source_thermal,                          &
        profile_list           = profile_list,                               &
        n_layer_stride         = nlayers+1,                                  &
        n_cloud_layer          = k,                                          &
        n_subcol_gen           = n_subcol_gen,                               &
        p_layer_1d             = pressure_in_wth(wth_1:wth_last),            &
        t_layer_1d             = temperature_in_wth(wth_1:wth_last),         &
        density_1d             = rho_in_wth(wth_1:wth_last),                 &
        cloud_frac_1d          = radiative_cloud_fraction(wth_1:wth_last),   &
        liq_frac_1d            = liquid_fraction(wth_1:wth_last),            &
        ice_frac_1d            = frozen_fraction(wth_1:wth_last),            &
        liq_mmr_1d             = mcl(wth_1:wth_last),                        &
        ice_mmr_1d             = mcf(wth_1:wth_last),                        &
        ice_nc_1d              = n_ice(wth_1:wth_last),                      &
        ice_conv_nc_1d         = conv_frozen_number(wth_1:wth_last),         &
        liq_dim_constant       = constant_droplet_effective_radius,          &
        liq_nc_1d              = cloud_drop_no_conc(wth_1:wth_last),         &
        conv_frac_1d           = radiative_conv_fraction(wth_1:wth_last),    &
        liq_conv_frac_1d       = conv_liquid_fraction(wth_1:wth_last),       &
        ice_conv_frac_1d       = conv_frozen_fraction(wth_1:wth_last),       &
        liq_conv_mmr_1d        = conv_liquid_mmr(wth_1:wth_last),            &
        ice_conv_mmr_1d        = conv_frozen_mmr(wth_1:wth_last),            &
        liq_conv_dim_constant  = constant_droplet_effective_radius,          &
        liq_conv_nc_1d         = cloud_drop_no_conc(wth_1:wth_last),         &
        liq_rsd_1d             = sigma_ml(wth_1:wth_last),                   &
        ice_rsd_1d             = sigma_mi(wth_1:wth_last),                   &
        cloud_vertical_decorr  = cloud_vertical_decorr,                      &
        conv_vertical_decorr   = cloud_vertical_decorr,                      &
        liq_dim_aparam         = liu_aparam,                                 &
        liq_dim_bparam         = liu_bparam,                                 &
        rand_seed              = rand_seed(twod_1:twod_last),                &
        l_mixing_ratio         = .true.,                                     &
        i_cloud_representation = i_cloud_representation,                     &
        i_overlap              = i_overlap,                                  &
        i_inhom                = ip_inhom_mcica,                             &
        i_drop_re              = i_drop_re,                                  &
        i_st_water             = i_cloud_liq_type_lw,                        &
        i_st_ice               = i_cloud_ice_type_lw,                        &
        i_cnv_water            = i_cloud_liq_type_lw,                        &
        i_cnv_ice              = i_cloud_ice_type_lw,                        &
        l_invert               = .true.,                                     &
        l_profile_last         = .true.)
    end if
  end do

  ! Set pointers for SW diagnostic output from runes call
  if (.not. associated(cloud_solar_extinction, empty_real_data)) then
    sw_diag%cloud_solar_extinction(0:nlayers, 1:n_profile) &
      => cloud_solar_extinction(wth_0:wth_last)
  else
    allocate( sw_diag%cloud_solar_extinction(0:nlayers, 1:n_profile) )
  end if

  ! Calculate cloud diagnostics from SW
  do k=0, nlayers
    profile_list = pack( [(l, l=1, n_profile)], &
                         n_cloud_layer(twod_1:twod_last) == k )
    n_profile_list = size(profile_list)
    if (n_profile_list > 0) then
      call runes(n_profile_list, nlayers, sw_diag,                           &
        spectrum_name          = 'sw',                                       &
        i_source               = ip_source_illuminate,                       &
        profile_list           = profile_list,                               &
        n_layer_stride         = nlayers+1,                                  &
        n_cloud_layer          = k,                                          &
        p_layer_1d             = pressure_in_wth(wth_1:wth_last),            &
        t_layer_1d             = temperature_in_wth(wth_1:wth_last),         &
        density_1d             = rho_in_wth(wth_1:wth_last),                 &
        cloud_frac_1d          = radiative_cloud_fraction(wth_1:wth_last),   &
        liq_frac_1d            = liquid_fraction(wth_1:wth_last),            &
        ice_frac_1d            = frozen_fraction(wth_1:wth_last),            &
        liq_mmr_1d             = mcl(wth_1:wth_last),                        &
        ice_mmr_1d             = mcf(wth_1:wth_last),                        &
        ice_nc_1d              = n_ice(wth_1:wth_last),                      &
        ice_conv_nc_1d         = conv_frozen_number(wth_1:wth_last),         &
        liq_dim_constant       = constant_droplet_effective_radius,          &
        liq_nc_1d              = cloud_drop_no_conc(wth_1:wth_last),         &
        conv_frac_1d           = radiative_conv_fraction(wth_1:wth_last),    &
        liq_conv_frac_1d       = conv_liquid_fraction(wth_1:wth_last),       &
        ice_conv_frac_1d       = conv_frozen_fraction(wth_1:wth_last),       &
        liq_conv_mmr_1d        = conv_liquid_mmr(wth_1:wth_last),            &
        ice_conv_mmr_1d        = conv_frozen_mmr(wth_1:wth_last),            &
        liq_conv_dim_constant  = constant_droplet_effective_radius,          &
        liq_conv_nc_1d         = cloud_drop_no_conc(wth_1:wth_last),         &
        cloud_vertical_decorr  = cloud_vertical_decorr,                      &
        conv_vertical_decorr   = cloud_vertical_decorr,                      &
        liq_dim_aparam         = liu_aparam,                                 &
        liq_dim_bparam         = liu_bparam,                                 &
        l_mixing_ratio         = .true.,                                     &
        i_cloud_representation = i_cloud_representation,                     &
        i_overlap              = i_overlap,                                  &
        i_inhom                = ip_inhom_homogeneous,                       &
        i_drop_re              = i_drop_re,                                  &
        i_st_water             = i_cloud_liq_type_sw,                        &
        i_st_ice               = i_cloud_ice_type_sw,                        &
        i_cnv_water            = i_cloud_liq_type_sw,                        &
        i_cnv_ice              = i_cloud_ice_type_sw,                        &
        l_invert               = .true.,                                     &
        l_profile_last         = .true.)
    end if
  end do

  ! Set COSP inputs
  p_full_levels(0:nlayers,1:n_profile) => pressure_in_wth(wth_0:wth_last)
  p_half_levels(1:nlayers,1:n_profile) => pressure_in_w3(w3_1:w3_last)
  hgt_full_levels(0:nlayers,1:n_profile) => height_wth(wth_0:wth_last)
  hgt_half_levels(1:nlayers,1:n_profile) => height_w3(w3_1:w3_last)
  density_full_levels(0:nlayers,1:n_profile) => rho_in_wth(wth_0:wth_last)
  t_full_levels(0:nlayers,1:n_profile) => temperature_in_wth(wth_0:wth_last)
  q_full_levels(0:nlayers,1:n_profile) => mv(wth_0:wth_last)
  cosp_crain(0:nlayers,1:n_profile) => conv_rain_3d(wth_0:wth_last)
  cosp_csnow(0:nlayers,1:n_profile) => conv_snow_3d(wth_0:wth_last)

  allocate(d_mass(nlayers, n_profile))
  allocate(t_surf(n_profile), p_surf(n_profile), hgt_surf(n_profile))
  if (.not. associated(sunlit_mask, empty_real_data)) then
    cosp_sunlit(1:n_profile) => sunlit_mask(twod_1:twod_last)
  else
    allocate(cosp_sunlit(n_profile))
  end if

  do l=1, n_profile
    d_mass(1, l) = density_full_levels(1, l) &
      * (hgt_half_levels(2, l) - hgt_full_levels(0, l))
    do k=2, nlayers-1
      d_mass(k, l) = density_full_levels(k, l) &
        * (hgt_half_levels(k+1, l) - hgt_half_levels(k, l))
    end do
    d_mass(nlayers, l) = density_full_levels(nlayers, l) * 2.0_r_def &
      * ( hgt_full_levels(nlayers, l) - hgt_half_levels(nlayers, l) )

    t_surf(l) = t_full_levels(0, l)
    p_surf(l) = p_full_levels(0, l)
    hgt_surf(l) = hgt_full_levels(0, l)

    if (lit_fraction(twod_1+l-1) > 0.0_r_def) then
      cosp_sunlit(l) = 1.0_r_def
    else
      cosp_sunlit(l) = 0.0_r_def
    end if
  end do

  ! Set COSP diagnostics
  if (.not. associated(calipso_low_cloud, empty_real_data)) then
    cosp_diag%cosp_calipso_low_level_cl(1:n_profile) &
      => calipso_low_cloud(twod_1:twod_last)
  end if
  if (.not. associated(calipso_low_cloud_mask, empty_real_data)) then
    cosp_diag%cosp_calipso_low_level_cl_mask(1:n_profile) &
      => calipso_low_cloud_mask(twod_1:twod_last)
  end if
  if (.not. associated(calipso_mid_cloud, empty_real_data)) then
    cosp_diag%cosp_calipso_mid_level_cl(1:n_profile) &
      => calipso_mid_cloud(twod_1:twod_last)
  end if
  if (.not. associated(calipso_mid_cloud_mask, empty_real_data)) then
    cosp_diag%cosp_calipso_mid_level_cl_mask(1:n_profile) &
      => calipso_mid_cloud_mask(twod_1:twod_last)
  end if
  if (.not. associated(calipso_high_cloud, empty_real_data)) then
    cosp_diag%cosp_calipso_high_level_cl(1:n_profile) &
      => calipso_high_cloud(twod_1:twod_last)
  end if
  if (.not. associated(calipso_high_cloud_mask, empty_real_data)) then
    cosp_diag%cosp_calipso_high_level_cl_mask(1:n_profile) &
      => calipso_high_cloud_mask(twod_1:twod_last)
  end if
  if (.not. associated(isccp_ctp_tau, empty_real_data)) then
    cosp_diag%cosp_ctp_tau_histogram(1:n_isccp_tau_bins, &
                                     1:n_isccp_pressure_bins, &
                                     1:n_profile) &
      => isccp_ctp_tau(ptau_1:ptau_last)
  end if
  if (.not. associated(calipso_cf_40_lvls_liq, empty_real_data)) then
    cosp_diag%cosp_calipso_cf_40_liq(1:n_cloudsat_levels, 1:n_profile) &
      => calipso_cf_40_lvls_liq(c40_1:c40_last)
  end if
  if (.not. associated(calipso_cf_40_lvls_ice, empty_real_data)) then
    cosp_diag%cosp_calipso_cf_40_ice(1:n_cloudsat_levels, 1:n_profile) &
      => calipso_cf_40_lvls_ice(c40_1:c40_last)
  end if
  if (.not. associated(calipso_cf_40_lvls_undet, empty_real_data)) then
    cosp_diag%cosp_calipso_cf_40_undet(1:n_cloudsat_levels, 1:n_profile) &
      => calipso_cf_40_lvls_undet(c40_1:c40_last)
  end if
  if (.not. associated(calipso_cf_40_lvls_mask, empty_real_data)) then
    cosp_diag%cosp_calipso_cf_40_mask(1:n_cloudsat_levels, 1:n_profile) &
      => calipso_cf_40_lvls_mask(c40_1:c40_last)
  end if
  if (.not. associated(calipso_cfad_sr_40_lvls, empty_real_data)) then
    cosp_diag%cosp_calipso_cfad_sr_40(1:n_backscatter_bins, &
                                      1:n_cloudsat_levels, &
                                      1:n_profile) &
      => calipso_cfad_sr_40_lvls(atb40_1:atb40_last)
  end if
  if (.not. associated(calipso_total_backscatter, empty_real_data)) then
    cosp_diag%cosp_calipso_tot_backscatter(0:nlayers, &
                                           1:n_subcol_gen, &
                                           1:n_profile) &
      => calipso_total_backscatter(subcol_0:subcol_last)
  end if

  ! Run COSP
  do k=0, nlayers
    profile_list = pack( [(l, l=1, n_profile)], &
                         n_cloud_layer(twod_1:twod_last) == k )
    n_profile_list = size(profile_list)
    if (k == 0) then
      n_subcol_cosp = 1
    else
      n_subcol_cosp = n_subcol_gen
    end if
    if (n_profile_list > 0) then
      call cosp(nlayers, n_profile_list, n_subcol_cosp, k, &
                lw_diag%n_subcol_cloud, &
                p_full_levels, &
                p_half_levels, &
                hgt_full_levels, &
                hgt_half_levels, &
                d_mass, t_full_levels, q_full_levels, &
                lw_diag%total_cloud_fraction, &
                lw_diag%liq_incloud_mmr, &
                lw_diag%ice_incloud_mmr, &
                cosp_crain, &
                cosp_csnow, &
                lw_diag%liq_dim, &
                lw_diag%ice_re, &
                sw_diag%cloud_solar_extinction, &
                lw_diag%cloud_thermal_absorptivity, &
                lw_diag%liq_part_frac, &
                lw_diag%ice_part_frac, &
                lw_diag%liq_subcol_scaling, &
                t_surf, p_surf, hgt_surf, &
                cosp_sunlit, &
                x1r, x1g, x2r, x2g, x4g, &
                cosp_diag, &
                profile_list = profile_list, &
                l_profile_last = .true. )
    end if
  end do

  if (associated(sunlit_mask, empty_real_data)) then
    deallocate(cosp_sunlit)
  end if
  deallocate(hgt_surf, p_surf, t_surf)
  deallocate(d_mass)
  if (associated(cloud_solar_extinction, empty_real_data)) then
    deallocate(sw_diag%cloud_solar_extinction)
  end if
  if (associated(cloud_thermal_absorptivity, empty_real_data)) then
    deallocate(lw_diag%cloud_thermal_absorptivity)
  end if
  deallocate(lw_diag%ice_incloud_mmr)
  deallocate(lw_diag%ice_part_frac)
  deallocate(lw_diag%ice_re)
  deallocate(lw_diag%liq_incloud_mmr)
  deallocate(lw_diag%liq_part_frac)
  deallocate(lw_diag%liq_dim)
  deallocate(lw_diag%total_cloud_fraction)
  deallocate(lw_diag%liq_subcol_scaling)
  deallocate(lw_diag%n_subcol_cloud)

end subroutine cosp_code
end module cosp_kernel_mod
