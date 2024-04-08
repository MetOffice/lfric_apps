!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to Socrates for shortwave fluxes (external illumination)

module sw_kernel_mod

use argument_mod,      only : arg_type, &
                              GH_FIELD, &
                              GH_REAL, GH_INTEGER, &
                              GH_READ, GH_WRITE, GH_READWRITE, &
                              DOMAIN, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5, &
                              ANY_DISCONTINUOUS_SPACE_6
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: sw_kernel_type
  private
  type(arg_type) :: meta_args(60) = (/ &
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    & ! sw_heating_rate_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_tile_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_blue_tile_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! rho_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! pressure_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! temperature_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! d_mass
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! layer_heat_capacity
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! stellar_irradiance_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! orographic_correction_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ozone
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mv
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mcl
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mci
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! n_ice
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_liquid_mmr
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_frozen_mmr
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! radiative_cloud_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! radiative_conv_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! frozen_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_frozen_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sigma_mc
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! cloud_drop_no_conc
    arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! rand_seed
    arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! n_cloud_layer
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_sw_direct_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_sw_diffuse_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sulphuric
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4), & ! aer_mix_ratio
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_sw_absorption
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_sw_scattering
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_sw_asymmetry
    ! Diagnostics (section radiation__)
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! cloud_top_re_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! cloud_top_weight_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! warm_cloud_top_re_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! warm_cloud_top_weight_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! cloud_extinction_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! cloud_weight_extinction_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! sw_aer_optical_depth_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_direct_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_down_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_up_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_direct_clear_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_down_clear_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_up_clear_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1)  & ! sw_up_clear_toa_rts
    /)
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: sw_code
end type

public :: sw_code

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                   Number of layers
!> @param[in]     n_profile                 Number of columns
!> @param[in,out] sw_heating_rate_rts       SW heating rate
!> @param[in,out] sw_down_surf_rts          SW downward surface flux
!> @param[in,out] sw_direct_surf_rts        SW unscattered surface flux
!> @param[in,out] sw_down_blue_surf_rts     SW blue downward surface flux
!> @param[in,out] sw_direct_blue_surf_rts   SW blue unscattered surface flux
!> @param[in,out] sw_up_surf_rts            SW upward surface flux
!> @param[in,out] sw_up_toa_rts             SW upward top-of-atmosphere flux
!> @param[in,out] sw_direct_toa_rts         SW unscattered top-of-atmosphere flux
!> @param[in,out] sw_up_tile_rts            SW upward tiled surface flux
!> @param[in,out] sw_up_blue_tile_rts       SW blue upward tiled surface flux
!> @param[in]     rho_in_wth                Density in wth space
!> @param[in]     pressure_in_wth           Pressure in wth space
!> @param[in]     temperature_in_wth        Temperature in wth space
!> @param[in]     d_mass                    Mass per square metre of radiation layers
!> @param[in]     layer_heat_capacity       Heat capacity of radiation layers
!> @param[in]     cos_zenith_angle_rts      Cosine of the stellar zenith angle
!> @param[in]     lit_fraction_rts          Lit fraction of the timestep
!> @param[in]     stellar_irradiance_rts    Stellar irradaince at the planet
!> @param[in]     orographic_correction_rts Orographic Correction
!> @param[in]     ozone                     Ozone field
!> @param[in]     mv                        Water vapour field
!> @param[in]     mcl                       Cloud liquid field
!> @param[in]     mci                       Cloud ice field
!> @param[in]     n_ice                     Ice number concentration
!> @param[in]     conv_liquid_mmr           Convective liquid gridbox MMR
!> @param[in]     conv_frozen_mmr           Convective frozen gridbox MMR
!> @param[in]     radiative_cloud_fraction  Large scale cloud fraction
!> @param[in]     radiative_conv_fraction   Convective cloud fraction
!> @param[in]     liquid_fraction           Liquid cloud fraction field
!> @param[in]     frozen_fraction           Frozen cloud fraction field
!> @param[in]     conv_liquid_fraction      Convective liquid cloud fraction
!> @param[in]     conv_frozen_fraction      Convective frozen cloud fraction
!> @param[in]     sigma_mc                  Fractional standard deviation of condensate
!> @param[in]     cloud_drop_no_conc        Cloud Droplet Number Concentration
!> @param[in]     rand_seed                 Random seed field for cloud generator
!> @param[in]     n_cloud_layer             Number of cloud layers
!> @param[in]     tile_fraction             Surface tile fractions
!> @param[in]     tile_sw_direct_albedo     SW direct tile albedos
!> @param[in]     tile_sw_diffuse_albedo    SW diffuse tile albedos
!> @param[in]     sulphuric                 Sulphuric acid aerosol
!> @param[in]     aer_mix_ratio             MODE aerosol mixing ratios
!> @param[in]     aer_sw_absorption         MODE aerosol SW absorption
!> @param[in]     aer_sw_scattering         MODE aerosol SW scattering
!> @param[in]     aer_sw_asymmetry          MODE aerosol SW asymmetry
!> @param[in,out] cloud_top_re_rts          Diagnostic: Cloud top effective radius
!> @param[in,out] cloud_top_weight_rts      Diagnostic: Cloud top weight
!> @param[in,out] warm_cloud_top_re_rts     Diagnostic: Warm cloud top effective radius
!> @param[in,out] warm_cloud_top_weight_rts Diagnostic: Warm cloud top weight
!> @param[in,out] cloud_extinction_rts      Diagnostic: Cloud extinction
!> @param[in,out] cloud_weight_extinction_rts     Diag: Weight for cloud extinction
!> @param[in,out] sw_aer_optical_depth_rts  Diagnostic: Aerosol optical depth in the visible
!> @param[in,out] sw_direct_rts             Diagnostic: SW direct flux on radiation levels
!> @param[in,out] sw_down_rts               Diagnostic: SW downwards flux on radiation levels
!> @param[in,out] sw_up_rts                 Diagnostic: SW upwards flux on radiation levels
!> @param[in,out] sw_direct_clear_rts       Diagnostic: Clear-sky SW direct flux on radiation levels
!> @param[in,out] sw_down_clear_rts         Diagnostic: Clear-sky SW downwards flux on radiation levels
!> @param[in,out] sw_up_clear_rts           Diagnostic: Clear-sky SW upwards flux on radiation levels
!> @param[in,out] sw_down_clear_surf_rts    Diagnostic: Clear-sky SW downwards surface flux
!> @param[in,out] sw_up_clear_surf_rts      Diagnostic: Clear-sky SW upwards surface flux
!> @param[in,out] sw_up_clear_toa_rts       Diagnostic: Clear-sky SW upwards top-of-atmosphere flux
!> @param[in]     ndf_wth                   No. DOFs per cell for wth space
!> @param[in]     undf_wth                  No. unique of DOFs for wth space
!> @param[in]     map_wth                   Dofmap for wth space column base cell
!> @param[in]     ndf_2d                    No. of DOFs per cell for 2D space
!> @param[in]     undf_2d                   No. unique of DOFs for 2D space
!> @param[in]     map_2d                    Dofmap for 2D space column base cell
!> @param[in]     ndf_tile                  Number of DOFs per cell for tiles
!> @param[in]     undf_tile                 Number of total DOFs for tiles
!> @param[in]     map_tile                  Dofmap for tile space column base cell
!> @param[in]     ndf_rtile                 No. of DOFs per cell for rtile space
!> @param[in]     undf_rtile                No. unique of DOFs for rtile space
!> @param[in]     map_rtile                 Dofmap for rtile space column base cell
!> @param[in]     ndf_mode                  No. of DOFs per cell for mode space
!> @param[in]     undf_mode                 No. unique of DOFs for mode space
!> @param[in]     map_mode                  Dofmap for mode space column base cell
!> @param[in]     ndf_rmode                 No. of DOFs per cell for rmode space
!> @param[in]     undf_rmode                No. unique of DOFs for rmode space
!> @param[in]     map_rmode                 Dofmap for rmode space column base cell
!> @param[in]     ndf_flux                  No. of DOFs per cell for flux space
!> @param[in]     undf_flux                 No. unique of DOFs for flux space
!> @param[in]     map_flux                  Dofmap for flux space column base cell
subroutine sw_code(nlayers, n_profile,                                         &
                   sw_heating_rate_rts, sw_down_surf_rts, sw_direct_surf_rts,  &
                   sw_down_blue_surf_rts, sw_direct_blue_surf_rts,             &
                   sw_up_surf_rts, sw_up_toa_rts, sw_direct_toa_rts,           &
                   sw_up_tile_rts, sw_up_blue_tile_rts,                        &
                   rho_in_wth, pressure_in_wth, temperature_in_wth,            &
                   d_mass, layer_heat_capacity,                                &
                   cos_zenith_angle_rts, lit_fraction_rts,                     &
                   stellar_irradiance_rts, orographic_correction_rts,          &
                   ozone, mv, mcl, mci, n_ice,                                 &
                   conv_liquid_mmr, conv_frozen_mmr,                           &
                   radiative_cloud_fraction, radiative_conv_fraction,          &
                   liquid_fraction, frozen_fraction,                           &
                   conv_liquid_fraction, conv_frozen_fraction,                 &
                   sigma_mc, cloud_drop_no_conc,                               &
                   rand_seed, n_cloud_layer,                                   &
                   tile_fraction,                                              &
                   tile_sw_direct_albedo, tile_sw_diffuse_albedo,              &
                   sulphuric, aer_mix_ratio,                                   &
                   aer_sw_absorption, aer_sw_scattering, aer_sw_asymmetry,     &
                   cloud_top_re_rts, cloud_top_weight_rts,                     &
                   warm_cloud_top_re_rts, warm_cloud_top_weight_rts,           &
                   cloud_extinction_rts, cloud_weight_extinction_rts,          &
                   sw_aer_optical_depth_rts,                                   &
                   sw_direct_rts, sw_down_rts, sw_up_rts,                      &
                   sw_direct_clear_rts, sw_down_clear_rts, sw_up_clear_rts,    &
                   sw_down_clear_surf_rts, sw_up_clear_surf_rts,               &
                   sw_up_clear_toa_rts,                                        &
                   ndf_wth, undf_wth, map_wth,                                 &
                   ndf_2d, undf_2d, map_2d,                                    &
                   ndf_tile, undf_tile, map_tile,                              &
                   ndf_rtile, undf_rtile, map_rtile,                           &
                   ndf_mode, undf_mode, map_mode,                              &
                   ndf_rmode, undf_rmode, map_rmode,                           &
                   ndf_flux, undf_flux, map_flux)

  use well_mixed_gases_config_mod, only: o2_mix_ratio
  use radiation_config_mod, only: &
    l_rayleigh_sw, &
    i_cloud_ice_type_sw, i_cloud_liq_type_sw, &
    cloud_vertical_decorr, &
    constant_droplet_effective_radius, liu_aparam, liu_bparam
  use aerosol_config_mod, only: l_radaer, sulphuric_strat_climatology, &
                                easyaerosol_sw
  use jules_control_init_mod, only: n_surf_tile
  use socrates_init_mod, only: n_sw_band, &
    i_cloud_representation, i_overlap, i_inhom, i_drop_re, l_orog
  use um_physics_init_mod, only: n_aer_mode_sw, mode_dimen, sw_band_mode
  use socrates_runes, only: runes, StrDiag, ip_source_illuminate
  use empty_data_mod, only: empty_real_data
  use gas_calc_all_mod, only: co2_mix_ratio_now,    &
                              n2o_mix_ratio_now,    &
                              ch4_mix_ratio_now

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, n_profile
  integer(i_def), intent(in) :: ndf_wth, undf_wth
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: ndf_rtile, undf_rtile
  integer(i_def), intent(in) :: ndf_mode, undf_mode
  integer(i_def), intent(in) :: ndf_rmode, undf_rmode
  integer(i_def), intent(in) :: ndf_flux, undf_flux

  integer(i_def), dimension(ndf_wth, n_profile),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d, n_profile),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile, n_profile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_rtile, n_profile), intent(in) :: map_rtile
  integer(i_def), dimension(ndf_mode, n_profile),  intent(in) :: map_mode
  integer(i_def), dimension(ndf_rmode, n_profile), intent(in) :: map_rmode
  integer(i_def), dimension(ndf_flux, n_profile),  intent(in) :: map_flux

  real(r_def), dimension(undf_wth),  intent(inout), target :: &
    sw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(inout), target :: &
    sw_down_surf_rts, sw_direct_surf_rts, &
    sw_up_surf_rts, sw_up_toa_rts, sw_direct_toa_rts
  real(r_def), dimension(undf_2d),   intent(inout), target :: &
    sw_down_blue_surf_rts, sw_direct_blue_surf_rts
  real(r_def), dimension(undf_tile), intent(inout), target :: &
    sw_up_tile_rts, sw_up_blue_tile_rts

  real(r_def), dimension(undf_wth), intent(in) :: &
    rho_in_wth, pressure_in_wth, temperature_in_wth, &
    d_mass, layer_heat_capacity, ozone, mv, mcl, mci, &
    n_ice, conv_liquid_mmr, conv_frozen_mmr, &
    radiative_cloud_fraction, radiative_conv_fraction, &
    liquid_fraction, frozen_fraction, &
    conv_liquid_fraction, conv_frozen_fraction, &
    sigma_mc, cloud_drop_no_conc

  integer(i_def), dimension(undf_2d), intent(in) :: rand_seed, n_cloud_layer

  real(r_def), dimension(undf_2d), intent(in) :: &
    cos_zenith_angle_rts, lit_fraction_rts, stellar_irradiance_rts, &
    orographic_correction_rts
  real(r_def), dimension(undf_tile),  intent(in) :: tile_fraction
  real(r_def), dimension(undf_rtile), intent(in) :: &
    tile_sw_direct_albedo, tile_sw_diffuse_albedo
  real(r_def), dimension(undf_wth),   intent(in) :: sulphuric
  real(r_def), dimension(undf_mode),  intent(in) :: aer_mix_ratio
  real(r_def), dimension(undf_rmode), intent(in) :: &
    aer_sw_absorption, aer_sw_scattering, aer_sw_asymmetry

  ! Conditional Diagnostics
  real(r_def), pointer, dimension(:), intent(inout) :: & ! 2d
    cloud_top_re_rts, cloud_top_weight_rts, &
    warm_cloud_top_re_rts, warm_cloud_top_weight_rts, &
    sw_down_clear_surf_rts, sw_up_clear_surf_rts, sw_up_clear_toa_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! wth
    cloud_extinction_rts, cloud_weight_extinction_rts, &
    sw_aer_optical_depth_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! flux
    sw_direct_rts, sw_down_rts, sw_up_rts, &
    sw_direct_clear_rts, sw_down_clear_rts, sw_up_clear_rts

  ! Local variables for the kernel
  integer(i_def) :: n_profile_list
  integer(i_def), allocatable :: unlit_list(:)
  integer(i_def), allocatable :: profile_list(:)
  integer(i_def) :: k, l
  integer(i_def) :: wth_0, wth_1, wth_last
  integer(i_def) :: tile_1, tile_last, rtile_1, rtile_last
  integer(i_def) :: mode_1, mode_last, rmode_1, rmode_last
  integer(i_def) :: flux_0, flux_last, twod_1, twod_last
  type(StrDiag)  :: sw_diag
  logical        :: l_aerosol_mode

  ! Set indexing
  wth_0 = map_wth(1,1)
  wth_1 = map_wth(1,1)+1
  wth_last = map_wth(1,1)+n_profile*(nlayers+1)-1
  tile_1 = map_tile(1,1)
  tile_last = map_tile(1,1)+n_profile*n_surf_tile-1
  rtile_1 = map_rtile(1,1)
  rtile_last = map_rtile(1,1)+n_profile*n_sw_band*n_surf_tile-1
  mode_1 = map_mode(1,1)+1
  mode_last = map_mode(1,1)+n_profile*(nlayers+1)*mode_dimen-1
  rmode_1 = map_rmode(1,1)+1
  rmode_last = map_rmode(1,1)+n_profile*(nlayers+1)*sw_band_mode-1
  flux_0 = map_flux(1,1)
  flux_last = map_flux(1,1)+n_profile*(nlayers+1)-1
  twod_1 = map_2d(1,1)
  twod_last = map_2d(1,1)+n_profile-1

  unlit_list = pack( [(l, l=1, n_profile)], &
                     lit_fraction_rts(twod_1:twod_last) <= 0.0_r_def )

  ! Set pointers for diagnostic output (using pointer bound remapping)
  ! Initialise fields to zero as unlit points will not be filled
  ! Always required:
  sw_diag%heating_rate(0:nlayers, 1:n_profile) &
                    => sw_heating_rate_rts(wth_0:wth_last)
  sw_diag%heating_rate(:, unlit_list) = 0.0_r_def
  sw_diag%flux_up_tile(1:n_surf_tile, 1:n_profile) &
                    => sw_up_tile_rts(tile_1:tile_last)
  sw_diag%flux_up_tile(:, unlit_list) = 0.0_r_def
  sw_diag%flux_up_blue_tile(1:n_surf_tile, 1:n_profile) &
                    => sw_up_blue_tile_rts(tile_1:tile_last)
  sw_diag%flux_up_blue_tile(:, unlit_list) = 0.0_r_def
  sw_diag%flux_direct_blue_surf(1:n_profile) &
                    => sw_direct_blue_surf_rts(twod_1:twod_last)
  sw_diag%flux_direct_blue_surf(unlit_list) = 0.0_r_def
  sw_diag%flux_down_blue_surf(1:n_profile) &
                    => sw_down_blue_surf_rts(twod_1:twod_last)
  sw_diag%flux_down_blue_surf(unlit_list) = 0.0_r_def
  sw_diag%flux_down_surf(1:n_profile) &
                    => sw_down_surf_rts(twod_1:twod_last)
  sw_diag%flux_down_surf(unlit_list) = 0.0_r_def
  sw_diag%flux_up_surf(1:n_profile) &
                    => sw_up_surf_rts(twod_1:twod_last)
  sw_diag%flux_up_surf(unlit_list) = 0.0_r_def
  sw_diag%flux_direct_surf(1:n_profile) &
                    => sw_direct_surf_rts(twod_1:twod_last)
  sw_diag%flux_direct_surf(unlit_list) = 0.0_r_def
  sw_diag%flux_up_toa(1:n_profile) &
                    => sw_up_toa_rts(twod_1:twod_last)
  sw_diag%flux_up_toa(unlit_list) = 0.0_r_def
  sw_diag%flux_direct_toa(1:n_profile) &
                    => sw_direct_toa_rts(twod_1:twod_last)
  sw_diag%flux_direct_toa(unlit_list) = 0.0_r_def

  ! Diagnosed on request:
  if (.not. associated(sw_direct_rts, empty_real_data)) then
    sw_diag%flux_direct(0:nlayers, 1:n_profile) &
                    => sw_direct_rts(flux_0:flux_last)
    sw_diag%flux_direct(:, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_down_rts, empty_real_data)) then
    sw_diag%flux_down(0:nlayers, 1:n_profile) &
                    => sw_down_rts(flux_0:flux_last)
    sw_diag%flux_down(:, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_rts, empty_real_data)) then
    sw_diag%flux_up(0:nlayers, 1:n_profile) &
                    => sw_up_rts(flux_0:flux_last)
    sw_diag%flux_up(:, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_direct_clear_rts, empty_real_data)) then
    sw_diag%flux_direct_clear(0:nlayers, 1:n_profile) &
                    => sw_direct_clear_rts(flux_0:flux_last)
    sw_diag%flux_direct_clear(:, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_down_clear_rts, empty_real_data)) then
    sw_diag%flux_down_clear(0:nlayers, 1:n_profile) &
                    => sw_down_clear_rts(flux_0:flux_last)
    sw_diag%flux_down_clear(:, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clear_rts, empty_real_data)) then
    sw_diag%flux_up_clear(0:nlayers, 1:n_profile) &
                    => sw_up_clear_rts(flux_0:flux_last)
    sw_diag%flux_up_clear(:, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_down_clear_surf_rts, empty_real_data)) then
    sw_diag%flux_down_clear_surf(1:n_profile) &
                    => sw_down_clear_surf_rts(twod_1:twod_last)
    sw_diag%flux_down_clear_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clear_surf_rts, empty_real_data)) then
    sw_diag%flux_up_clear_surf(1:n_profile) &
                    => sw_up_clear_surf_rts(twod_1:twod_last)
    sw_diag%flux_up_clear_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clear_toa_rts, empty_real_data)) then
    sw_diag%flux_up_clear_toa(1:n_profile) &
                    => sw_up_clear_toa_rts(twod_1:twod_last)
    sw_diag%flux_up_clear_toa(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(cloud_top_re_rts, empty_real_data) .and. &
      .not. associated(cloud_top_weight_rts, empty_real_data)) then
    sw_diag%cloud_top_liq_dim(1:n_profile) &
                    => cloud_top_re_rts(twod_1:twod_last)
    sw_diag%cloud_top_liq_dim(unlit_list) = 0.0_r_def
    sw_diag%cloud_top_liq_weight(1:n_profile) &
                    => cloud_top_weight_rts(twod_1:twod_last)
    sw_diag%cloud_top_liq_weight(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(warm_cloud_top_re_rts, empty_real_data) .and. &
      .not. associated(warm_cloud_top_weight_rts, empty_real_data)) then
    sw_diag%cloud_top_warm_liq_dim(1:n_profile) &
                    => warm_cloud_top_re_rts(twod_1:twod_last)
    sw_diag%cloud_top_warm_liq_dim(unlit_list) = 0.0_r_def
    sw_diag%cloud_top_warm_liq_weight(1:n_profile) &
                    => warm_cloud_top_weight_rts(twod_1:twod_last)
    sw_diag%cloud_top_warm_liq_weight(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(cloud_extinction_rts, empty_real_data)) then
    sw_diag%cloud_extinction(0:nlayers, 1:n_profile) &
                    => cloud_extinction_rts(wth_0:wth_last)
    sw_diag%cloud_extinction(:, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(cloud_weight_extinction_rts, empty_real_data)) then
    sw_diag%cloud_weight_extinction(0:nlayers, 1:n_profile) &
                    => cloud_weight_extinction_rts(wth_0:wth_last)
    sw_diag%cloud_weight_extinction(:, unlit_list) = 0.0_r_def
  end if

  ! Aerosol optical depth for SW band 3 is output (505-690nm). Once the
  ! diagnostic infrastructure is in place the band(s) may be user defined.
  if (.not. associated(sw_aer_optical_depth_rts, empty_real_data)) then
    sw_diag%aerosol_optical_depth(0:nlayers, 3:3, 1:n_profile) &
                    => sw_aer_optical_depth_rts(wth_0:wth_last)
    sw_diag%aerosol_optical_depth(:, :, unlit_list) = 0.0_r_def
  end if

  ! If radaer or easyaerosol are running, socrates includes aerosol modes
  if (l_radaer .or. easyaerosol_sw) then
    l_aerosol_mode = .true.
  else
    l_aerosol_mode = .false.
  end if

  do k=0, nlayers
    profile_list = pack( [(l, l=1, n_profile)], &
                         lit_fraction_rts(twod_1:twod_last) > 0.0_r_def &
                         .and. n_cloud_layer(twod_1:twod_last) == k )
    n_profile_list = size(profile_list)
    if (n_profile_list > 0) then
      ! Calculate the SW fluxes (RUN the Edwards-Slingo two-stream solver)
      call runes(n_profile_list, nlayers, sw_diag,                           &
        spectrum_name          = 'sw',                                       &
        i_source               = ip_source_illuminate,                       &
        profile_list           = profile_list,                               &
        n_layer_stride         = nlayers+1,                                  &
        n_cloud_layer          = k,                                          &
        p_layer_1d             = pressure_in_wth(wth_1:wth_last),            &
        t_layer_1d             = temperature_in_wth(wth_1:wth_last),         &
        mass_1d                = d_mass(wth_1:wth_last),                     &
        density_1d             = rho_in_wth(wth_1:wth_last),                 &
        h2o_1d                 = mv(wth_1:wth_last),                         &
        o3_1d                  = ozone(wth_1:wth_last),                      &
        co2_mix_ratio          = co2_mix_ratio_now,                          &
        n2o_mix_ratio          = n2o_mix_ratio_now,                          &
        ch4_mix_ratio          = ch4_mix_ratio_now,                          &
        o2_mix_ratio           = o2_mix_ratio,                               &
        cos_zenith_angle       = cos_zenith_angle_rts(twod_1:twod_last),     &
        solar_irrad            = stellar_irradiance_rts(twod_1:twod_last),   &
        l_orog                 = l_orog,                                     &
        orog_corr              = orographic_correction_rts(twod_1:twod_last),&
        n_tile                 = n_surf_tile,                                &
        frac_tile_1d           = tile_fraction(tile_1:tile_last),            &
        albedo_diff_tile_1d    = tile_sw_diffuse_albedo(rtile_1:rtile_last), &
        albedo_dir_tile_1d     = tile_sw_direct_albedo(rtile_1:rtile_last),  &
        cloud_frac_1d          = radiative_cloud_fraction(wth_1:wth_last),   &
        liq_frac_1d            = liquid_fraction(wth_1:wth_last),            &
        ice_frac_1d            = frozen_fraction(wth_1:wth_last),            &
        liq_mmr_1d             = mcl(wth_1:wth_last),                        &
        ice_mmr_1d             = mci(wth_1:wth_last),                        &
        ice_nc_1d              = n_ice(wth_1:wth_last),                      &
        ice_conv_nc_1d         = n_ice(wth_1:wth_last),                      &
        liq_dim_constant       = constant_droplet_effective_radius,          &
        liq_nc_1d              = cloud_drop_no_conc(wth_1:wth_last),         &
        conv_frac_1d           = radiative_conv_fraction(wth_1:wth_last),    &
        liq_conv_frac_1d       = conv_liquid_fraction(wth_1:wth_last),       &
        ice_conv_frac_1d       = conv_frozen_fraction(wth_1:wth_last),       &
        liq_conv_mmr_1d        = conv_liquid_mmr(wth_1:wth_last),            &
        ice_conv_mmr_1d        = conv_frozen_mmr(wth_1:wth_last),            &
        liq_conv_dim_constant  = constant_droplet_effective_radius,          &
        liq_conv_nc_1d         = cloud_drop_no_conc(wth_1:wth_last),         &
        liq_rsd_1d             = sigma_mc(wth_1:wth_last),                   &
        ice_rsd_1d             = sigma_mc(wth_1:wth_last),                   &
        cloud_vertical_decorr  = cloud_vertical_decorr,                      &
        conv_vertical_decorr   = cloud_vertical_decorr,                      &
        liq_dim_aparam         = liu_aparam,                                 &
        liq_dim_bparam         = liu_bparam,                                 &
        rand_seed              = rand_seed(twod_1:twod_last),                &
        layer_heat_capacity_1d = layer_heat_capacity(wth_1:wth_last),        &
        l_rayleigh             = l_rayleigh_sw,                              &
        l_mixing_ratio         = .true.,                                     &
        i_cloud_representation = i_cloud_representation,                     &
        i_overlap              = i_overlap,                                  &
        i_inhom                = i_inhom,                                    &
        i_drop_re              = i_drop_re,                                  &
        i_st_water             = i_cloud_liq_type_sw,                        &
        i_st_ice               = i_cloud_ice_type_sw,                        &
        i_cnv_water            = i_cloud_liq_type_sw,                        &
        i_cnv_ice              = i_cloud_ice_type_sw,                        &
        l_sulphuric            = sulphuric_strat_climatology,                &
        sulphuric_1d           = sulphuric(wth_1:wth_last),                  &
        l_aerosol_mode         = l_aerosol_mode,                             &
        n_aer_mode             = n_aer_mode_sw,                              &
        aer_mix_ratio_1d       = aer_mix_ratio(mode_1:mode_last),            &
        aer_absorption_1d      = aer_sw_absorption(rmode_1:rmode_last),      &
        aer_scattering_1d      = aer_sw_scattering(rmode_1:rmode_last),      &
        aer_asymmetry_1d       = aer_sw_asymmetry(rmode_1:rmode_last),       &
        l_invert               = .true.,                                     &
        l_profile_last         = .true.)
    end if
  end do

end subroutine sw_code
end module sw_kernel_mod
