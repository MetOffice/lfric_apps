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
                              ANY_DISCONTINUOUS_SPACE_6, &
                              ANY_DISCONTINUOUS_SPACE_7, &
                              ANY_DISCONTINUOUS_SPACE_8
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
  type(arg_type) :: meta_args(99) = (/ &
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
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! h2o
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! co2
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! o3
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! n2o
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! co
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ch4
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! o2
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! so2
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! nh3
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! n2
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! h2
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! he
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! hcn
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! cs
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! k
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! li
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! na
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! rb
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! tio
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! vo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mcl
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mcf
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! n_ice
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_frozen_number
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_liquid_mmr
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_frozen_mmr
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! radiative_cloud_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! radiative_conv_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! frozen_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! conv_frozen_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sigma_ml
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sigma_mi
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
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7), & ! sw_direct_band_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7), & ! sw_down_band_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7), & ! sw_up_band_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7), & ! sw_direct_clear_band_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7), & ! sw_down_clear_band_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7), & ! sw_up_clear_band_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clear_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_clean_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clean_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clean_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_clear_clean_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clear_clean_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clear_clean_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_uv_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_uv_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_uv_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_uv_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_uv_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_uv_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_8)  & ! photolysis_rates_rts
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
!> @param[in]     h2o                       Water vapour
!> @param[in]     co2                       Carbon dioxide
!> @param[in]     o3                        Ozone
!> @param[in]     n2o                       Dinitrogen oxide
!> @param[in]     co                        Carbon monoxide
!> @param[in]     ch4                       Methane
!> @param[in]     o2                        Oxygen
!> @param[in]     so2                       Sulphur dioxide
!> @param[in]     nh3                       Ammonia
!> @param[in]     n2                        Nitrogen
!> @param[in]     h2                        Hydrogen
!> @param[in]     he                        Helium
!> @param[in]     hcn                       Hydrogen cyanide
!> @param[in]     cs                        Cesium
!> @param[in]     potassium                 Potassium
!> @param[in]     li                        Lithium
!> @param[in]     na                        Sodium
!> @param[in]     rb                        Rubidium
!> @param[in]     tio                       Titanium oxide
!> @param[in]     vo                        Vanadium oxide
!> @param[in]     mcl                       Cloud liquid field
!> @param[in]     mcf                       Cloud ice field
!> @param[in]     n_ice                     Ice number concentration
!> @param[in]     conv_frozen_number        Convective Ice number concentration
!> @param[in]     conv_liquid_mmr           Convective liquid gridbox MMR
!> @param[in]     conv_frozen_mmr           Convective frozen gridbox MMR
!> @param[in]     radiative_cloud_fraction  Large scale cloud fraction
!> @param[in]     radiative_conv_fraction   Convective cloud fraction
!> @param[in]     liquid_fraction           Liquid cloud fraction field
!> @param[in]     frozen_fraction           Frozen cloud fraction field
!> @param[in]     conv_liquid_fraction      Convective liquid cloud fraction
!> @param[in]     conv_frozen_fraction      Convective frozen cloud fraction
!> @param[in]     sigma_ml                  Fractional standard deviation of liquid condensate
!> @param[in]     sigma_mi                  Fractional standard deviation of ice condensate
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
!> @param[in,out] sw_direct_band_rts        Diagnostic: SW direct flux on radiation levels and bands
!> @param[in,out] sw_down_band_rts          Diagnostic: SW downwards flux on radiation levels and bands
!> @param[in,out] sw_up_band_rts            Diagnostic: SW upwards flux on radiation levels and bands
!> @param[in,out] sw_direct_clear_band_rts  Diagnostic: Clear-sky SW direct flux on radiation levels and bands
!> @param[in,out] sw_down_clear_band_rts    Diagnostic: Clear-sky SW down flux on radiation levels and bands
!> @param[in,out] sw_up_clear_band_rts      Diagnostic: Clear-sky SW up flux on radiation levels and bands
!> @param[in,out] sw_down_clear_surf_rts    Diagnostic: Clear-sky SW downwards surface flux
!> @param[in,out] sw_up_clear_surf_rts      Diagnostic: Clear-sky SW upwards surface flux
!> @param[in,out] sw_up_clear_toa_rts       Diagnostic: Clear-sky SW upwards top-of-atmosphere flux
!> @param[in,out] sw_down_clean_surf_rts    Diagnostic: Clean-air SW downwards surface flux
!> @param[in,out] sw_up_clean_surf_rts      Diagnostic: Clean-air SW upwards surface flux
!> @param[in,out] sw_up_clean_toa_rts       Diagnostic: Clean-air SW upwards top-of-atmosphere flux
!> @param[in,out] sw_down_clear_clean_surf_rts    Diag: Clear-clean SW downwards surface flux
!> @param[in,out] sw_up_clear_clean_surf_rts      Diag: Clear-clean SW upwards surface flux
!> @param[in,out] sw_up_clear_clean_toa_rts       Diag: Clear-clean SW upwards top-of-atmosphere flux
!> @param[in,out] sw_down_uv_surf_rts       Diagnostic: UV downwards surface flux
!> @param[in,out] sw_down_uv_clear_surf_rts Diagnostic: Clear-sky UV downwards surface flux
!> @param[in,out] sw_direct_uv_surf_rts     Diagnostic: UV direct surface flux
!> @param[in,out] sw_direct_uv_clear_surf_rts     Diag: Clear-sky UV direct surface flux
!> @param[in,out] sw_up_uv_surf_rts         Diagnostic: UV upwards surface flux
!> @param[in,out] sw_up_uv_clear_surf_rts   Diagnostic: Clear-sky UV upwards surface flux
!> @param[in,out] photolysis_rates_rts      Diagnostic: Photolysis rates
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
!> @param[in]     ndf_bflux                 No. of DOFs per cell for bflux space
!> @param[in]     undf_bflux                No. unique of DOFs for bflux space
!> @param[in]     map_bflux                 Dofmap for bflux space column base cell
!> @param[in]     ndf_ph                    No. of DOFs per cell for ph space
!> @param[in]     undf_ph                   No. unique of DOFs for ph space
!> @param[in]     map_ph                    Dofmap for ph space column base cell
subroutine sw_code(nlayers, n_profile, &
    sw_heating_rate_rts, sw_down_surf_rts, sw_direct_surf_rts, &
    sw_down_blue_surf_rts, sw_direct_blue_surf_rts, &
    sw_up_surf_rts, sw_up_toa_rts, sw_direct_toa_rts, &
    sw_up_tile_rts, sw_up_blue_tile_rts, &
    rho_in_wth, pressure_in_wth, temperature_in_wth, &
    d_mass, layer_heat_capacity, &
    cos_zenith_angle_rts, lit_fraction_rts, &
    stellar_irradiance_rts, orographic_correction_rts, &
    h2o, co2, o3, n2o, co, ch4, o2, so2, nh3, n2, h2, he, hcn, &
    cs, potassium, li, na, rb, tio, vo, &
    mcl, mcf, n_ice, conv_frozen_number, &
    conv_liquid_mmr, conv_frozen_mmr, &
    radiative_cloud_fraction, radiative_conv_fraction, &
    liquid_fraction, frozen_fraction, &
    conv_liquid_fraction, conv_frozen_fraction, &
    sigma_ml, sigma_mi, cloud_drop_no_conc, &
    rand_seed, n_cloud_layer, &
    tile_fraction, &
    tile_sw_direct_albedo, tile_sw_diffuse_albedo, &
    sulphuric, aer_mix_ratio, &
    aer_sw_absorption, aer_sw_scattering, aer_sw_asymmetry, &
    ! Conditional diagnostics
    cloud_top_re_rts, cloud_top_weight_rts, &
    warm_cloud_top_re_rts, warm_cloud_top_weight_rts, &
    cloud_extinction_rts, cloud_weight_extinction_rts, &
    sw_aer_optical_depth_rts, &
    sw_direct_rts, sw_down_rts, sw_up_rts, &
    sw_direct_clear_rts, sw_down_clear_rts, sw_up_clear_rts, &
    sw_direct_band_rts, sw_down_band_rts, sw_up_band_rts, &
    sw_direct_clear_band_rts, sw_down_clear_band_rts, sw_up_clear_band_rts, &
    sw_down_clear_surf_rts, sw_up_clear_surf_rts, sw_up_clear_toa_rts, &
    sw_down_clean_surf_rts, sw_up_clean_surf_rts, sw_up_clean_toa_rts, &
    sw_down_clear_clean_surf_rts, sw_up_clear_clean_surf_rts, &
    sw_up_clear_clean_toa_rts, &
    sw_down_uv_surf_rts, sw_down_uv_clear_surf_rts, &
    sw_direct_uv_surf_rts, sw_direct_uv_clear_surf_rts, &
    sw_up_uv_surf_rts, sw_up_uv_clear_surf_rts, &
    photolysis_rates_rts, &
    ndf_wth, undf_wth, map_wth, &
    ndf_2d, undf_2d, map_2d, &
    ndf_tile, undf_tile, map_tile, &
    ndf_rtile, undf_rtile, map_rtile, &
    ndf_mode, undf_mode, map_mode, &
    ndf_rmode, undf_rmode, map_rmode, &
    ndf_flux, undf_flux, map_flux, &
    ndf_bflux, undf_bflux, map_bflux, &
    ndf_ph, undf_ph, map_ph)

  use radiation_config_mod, only: &
    l_rayleigh_sw, &
    i_cloud_ice_type_sw, i_cloud_liq_type_sw, &
    cloud_vertical_decorr, &
    constant_droplet_effective_radius, liu_aparam, liu_bparam
  use aerosol_config_mod, only: l_radaer, sulphuric_strat_climatology, &
                                easyaerosol_sw
  use jules_control_init_mod, only: n_surf_tile
  use socrates_init_mod, only: n_sw_band, n_pathway, sw_550nm, &
    i_cloud_representation, i_overlap, i_inhom, i_cloud_entrapment, &
    i_drop_re, l_orog
  use um_physics_init_mod, only: n_aer_mode_sw, mode_dimen, sw_band_mode
  use socrates_runes, only: runes, StrDiag, ip_source_illuminate
  use empty_data_mod, only: empty_real_data
  use gas_calc_all_mod, only: &
    cfc11_mix_ratio_now,  &
    cfc113_mix_ratio_now, &
    cfc12_mix_ratio_now,  &
    ch4_mix_ratio_now, ch4_well_mixed, &
    co_mix_ratio_now, co_well_mixed, &
    co2_mix_ratio_now, co2_well_mixed, &
    cs_mix_ratio_now, cs_well_mixed, &
    h2_mix_ratio_now, h2_well_mixed, &
    h2o_mix_ratio_now, h2o_well_mixed, &
    hcfc22_mix_ratio_now, &
    hcn_mix_ratio_now, hcn_well_mixed, &
    he_mix_ratio_now, he_well_mixed, &
    hfc134a_mix_ratio_now, &
    k_mix_ratio_now, k_well_mixed, &
    li_mix_ratio_now, li_well_mixed, &
    n2_mix_ratio_now, n2_well_mixed, &
    n2o_mix_ratio_now, n2o_well_mixed, &
    na_mix_ratio_now, na_well_mixed, &
    nh3_mix_ratio_now, nh3_well_mixed, &
    o2_mix_ratio_now, o2_well_mixed, &
    o3_mix_ratio_now, o3_well_mixed, &
    rb_mix_ratio_now, rb_well_mixed, &
    so2_mix_ratio_now, so2_well_mixed, &
    tio_mix_ratio_now, tio_well_mixed, &
    vo_mix_ratio_now, vo_well_mixed

  !$ use omp_lib, only: omp_get_max_threads

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
  integer(i_def), intent(in) :: ndf_bflux, undf_bflux
  integer(i_def), intent(in) :: ndf_ph, undf_ph

  integer(i_def), dimension(ndf_wth, n_profile),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d, n_profile),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile, n_profile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_rtile, n_profile), intent(in) :: map_rtile
  integer(i_def), dimension(ndf_mode, n_profile),  intent(in) :: map_mode
  integer(i_def), dimension(ndf_rmode, n_profile), intent(in) :: map_rmode
  integer(i_def), dimension(ndf_flux, n_profile),  intent(in) :: map_flux
  integer(i_def), dimension(ndf_bflux, n_profile), intent(in) :: map_bflux
  integer(i_def), dimension(ndf_ph, n_profile),    intent(in) :: map_ph

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
    d_mass, layer_heat_capacity, mcl, mcf, &
    n_ice, conv_frozen_number, conv_liquid_mmr, conv_frozen_mmr, &
    radiative_cloud_fraction, radiative_conv_fraction, &
    liquid_fraction, frozen_fraction, &
    conv_liquid_fraction, conv_frozen_fraction, &
    sigma_ml, sigma_mi, cloud_drop_no_conc, &
    h2o, co2, o3, n2o, co, ch4, o2, so2, nh3, n2, h2, he, hcn, &
    cs, potassium, li, na, rb, tio, vo

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
    sw_down_clear_surf_rts, sw_up_clear_surf_rts, sw_up_clear_toa_rts, &
    sw_down_clean_surf_rts, sw_up_clean_surf_rts, sw_up_clean_toa_rts, &
    sw_down_clear_clean_surf_rts, sw_up_clear_clean_surf_rts, &
    sw_up_clear_clean_toa_rts, &
    sw_down_uv_surf_rts, sw_down_uv_clear_surf_rts, &
    sw_direct_uv_surf_rts, sw_direct_uv_clear_surf_rts, &
    sw_up_uv_surf_rts, sw_up_uv_clear_surf_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! wth
    cloud_extinction_rts, cloud_weight_extinction_rts, &
    sw_aer_optical_depth_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! flux
    sw_direct_rts, sw_down_rts, sw_up_rts, &
    sw_direct_clear_rts, sw_down_clear_rts, sw_up_clear_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! bflux
    sw_direct_band_rts, sw_down_band_rts, sw_up_band_rts, &
    sw_direct_clear_band_rts, sw_down_clear_band_rts, sw_up_clear_band_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! ph
    photolysis_rates_rts

  ! Local variables for the kernel
  integer(i_def) :: n_profile_list, n_profile_list_seg
  integer(i_def), allocatable :: unlit_list(:)
  integer(i_def), allocatable :: profile_list(:)
  integer(i_def) :: k, l, ll
  integer(i_def) :: wth_0, wth_1, wth_last
  integer(i_def) :: tile_1, tile_last, rtile_1, rtile_last
  integer(i_def) :: mode_1, mode_last, rmode_1, rmode_last
  integer(i_def) :: flux_0, flux_last, bflux_0, bflux_last
  integer(i_def) :: ph_0, ph_last, twod_1, twod_last
  type(StrDiag)  :: sw_diag
  logical        :: l_aerosol_mode

  ! Segmentation variables for threading call to Socrates
  integer(i_def) :: max_threads, soc_sw_block, seg_start, seg_end

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
  bflux_0 = map_bflux(1,1)
  bflux_last = map_bflux(1,1)+n_profile*n_sw_band*(nlayers+1)-1
  ph_0 = map_ph(1,1)
  ph_last = map_ph(1,1)+n_profile*n_pathway*(nlayers+1)-1
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

  if (.not. associated(sw_direct_band_rts, empty_real_data)) then
    sw_diag%flux_direct_band(0:nlayers, 1:n_sw_band, 1:n_profile) &
                    => sw_direct_band_rts(bflux_0:bflux_last)
    sw_diag%flux_direct_band(:, :, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_down_band_rts, empty_real_data)) then
    sw_diag%flux_down_band(0:nlayers, 1:n_sw_band, 1:n_profile) &
                    => sw_down_band_rts(bflux_0:bflux_last)
    sw_diag%flux_down_band(:, :, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_band_rts, empty_real_data)) then
    sw_diag%flux_up_band(0:nlayers, 1:n_sw_band, 1:n_profile) &
                    => sw_up_band_rts(bflux_0:bflux_last)
    sw_diag%flux_up_band(:, :, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_direct_clear_band_rts, empty_real_data)) then
    sw_diag%flux_direct_clear_band(0:nlayers, 1:n_sw_band, 1:n_profile) &
                    => sw_direct_clear_band_rts(bflux_0:bflux_last)
    sw_diag%flux_direct_clear_band(:, :, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_down_clear_band_rts, empty_real_data)) then
    sw_diag%flux_down_clear_band(0:nlayers, 1:n_sw_band, 1:n_profile) &
                    => sw_down_clear_band_rts(bflux_0:bflux_last)
    sw_diag%flux_down_clear_band(:, :, unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clear_band_rts, empty_real_data)) then
    sw_diag%flux_up_clear_band(0:nlayers, 1:n_sw_band, 1:n_profile) &
                    => sw_up_clear_band_rts(bflux_0:bflux_last)
    sw_diag%flux_up_clear_band(:, :, unlit_list) = 0.0_r_def
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

  if (.not. associated(sw_down_clean_surf_rts, empty_real_data)) then
    sw_diag%flux_down_clean_surf(1:n_profile) &
                    => sw_down_clean_surf_rts(twod_1:twod_last)
    sw_diag%flux_down_clean_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clean_surf_rts, empty_real_data)) then
    sw_diag%flux_up_clean_surf(1:n_profile) &
                    => sw_up_clean_surf_rts(twod_1:twod_last)
    sw_diag%flux_up_clean_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clean_toa_rts, empty_real_data)) then
    sw_diag%flux_up_clean_toa(1:n_profile) &
                    => sw_up_clean_toa_rts(twod_1:twod_last)
    sw_diag%flux_up_clean_toa(unlit_list) = 0.0_r_def
  end if

  if (.not. associated(sw_down_clear_clean_surf_rts, empty_real_data)) then
    sw_diag%flux_down_clear_clean_surf(1:n_profile) &
                    => sw_down_clear_clean_surf_rts(twod_1:twod_last)
    sw_diag%flux_down_clear_clean_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clear_clean_surf_rts, empty_real_data)) then
    sw_diag%flux_up_clear_clean_surf(1:n_profile) &
                    => sw_up_clear_clean_surf_rts(twod_1:twod_last)
    sw_diag%flux_up_clear_clean_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_clear_clean_toa_rts, empty_real_data)) then
    sw_diag%flux_up_clear_clean_toa(1:n_profile) &
                    => sw_up_clear_clean_toa_rts(twod_1:twod_last)
    sw_diag%flux_up_clear_clean_toa(unlit_list) = 0.0_r_def
  end if

  if (.not. associated(sw_down_uv_surf_rts, empty_real_data)) then
    sw_diag%flux_down_uv_surf(1:n_profile) &
                    => sw_down_uv_surf_rts(twod_1:twod_last)
    sw_diag%flux_down_uv_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_down_uv_clear_surf_rts, empty_real_data)) then
    sw_diag%flux_down_uv_clear_surf(1:n_profile) &
                    => sw_down_uv_clear_surf_rts(twod_1:twod_last)
    sw_diag%flux_down_uv_clear_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_direct_uv_surf_rts, empty_real_data)) then
    sw_diag%flux_direct_uv_surf(1:n_profile) &
                    => sw_direct_uv_surf_rts(twod_1:twod_last)
    sw_diag%flux_direct_uv_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_direct_uv_clear_surf_rts, empty_real_data)) then
    sw_diag%flux_direct_uv_clear_surf(1:n_profile) &
                    => sw_direct_uv_clear_surf_rts(twod_1:twod_last)
    sw_diag%flux_direct_uv_clear_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_uv_surf_rts, empty_real_data)) then
    sw_diag%flux_up_uv_surf(1:n_profile) &
                    => sw_up_uv_surf_rts(twod_1:twod_last)
    sw_diag%flux_up_uv_surf(unlit_list) = 0.0_r_def
  end if
  if (.not. associated(sw_up_uv_clear_surf_rts, empty_real_data)) then
    sw_diag%flux_up_uv_clear_surf(1:n_profile) &
                    => sw_up_uv_clear_surf_rts(twod_1:twod_last)
    sw_diag%flux_up_uv_clear_surf(unlit_list) = 0.0_r_def
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

  if (.not. associated(photolysis_rates_rts, empty_real_data)) then
    sw_diag%photolysis_rate(0:nlayers, 1:n_pathway, 1:n_profile) &
                    => photolysis_rates_rts(ph_0:ph_last)
    sw_diag%photolysis_rate(:, :, unlit_list) = 0.0_r_def
  end if

  ! Aerosol optical depth for SW band at 550nm
  if (.not. associated(sw_aer_optical_depth_rts, empty_real_data)) then
    sw_diag%aerosol_optical_depth(0:nlayers, sw_550nm:sw_550nm, 1:n_profile) &
                    => sw_aer_optical_depth_rts(wth_0:wth_last)
    sw_diag%aerosol_optical_depth(:, :, unlit_list) = 0.0_r_def
  end if

  ! If radaer or easyaerosol are running, socrates includes aerosol modes
  if (l_radaer .or. easyaerosol_sw) then
    l_aerosol_mode = .true.
  else
    l_aerosol_mode = .false.
  end if

  do k = 0, nlayers
    profile_list = pack( [(l, l=1, n_profile)], &
                         lit_fraction_rts(twod_1:twod_last) > 0.0_r_def &
                         .and. n_cloud_layer(twod_1:twod_last) == k )
    n_profile_list = size(profile_list)
    if (n_profile_list > 0) then

      max_threads = 1
      !$ max_threads = omp_get_max_threads()
      soc_sw_block = ceiling(real(n_profile_list)/max_threads)

      ! The number of indices to run through per k step can become quite small
      ! this check ensures that should the generated block exceed the length
      ! of n_profile_list it is set to n_profile_list
      soc_sw_block = min(soc_sw_block, n_profile_list)

      !$OMP PARALLEL DEFAULT(SHARED)                                           &
      !$OMP PRIVATE(ll, seg_start, seg_end, n_profile_list_seg)
      !$OMP do SCHEDULE(STATIC)
      do ll=1, n_profile_list, soc_sw_block
        seg_start = ll

        if ((seg_start + soc_sw_block) -1 >= n_profile_list) then
          seg_end = n_profile_list
          n_profile_list_seg = (seg_end - seg_start) + 1
        else
          seg_end = (seg_start + soc_sw_block) - 1
          n_profile_list_seg = soc_sw_block
        end if

        ! Calculate the SW fluxes (RUN the Edwards-Slingo two-stream solver)
        call runes(n_profile_list_seg, nlayers, sw_diag,                       &
          spectrum_name          = 'sw',                                       &
          i_source               = ip_source_illuminate,                       &
          profile_list           = profile_list(seg_start:seg_end),            &
          n_layer_stride         = nlayers+1,                                  &
          n_cloud_layer          = k,                                          &
          p_layer_1d             = pressure_in_wth(wth_1:wth_last),            &
          t_layer_1d             = temperature_in_wth(wth_1:wth_last),         &
          mass_1d                = d_mass(wth_1:wth_last),                     &
          density_1d             = rho_in_wth(wth_1:wth_last),                 &
          ch4_1d                 = ch4(wth_1:wth_last),                        &
          co_1d                  = co(wth_1:wth_last),                         &
          co2_1d                 = co2(wth_1:wth_last),                        &
          cs_1d                  = cs(wth_1:wth_last),                         &
          h2_1d                  = h2(wth_1:wth_last),                         &
          h2o_1d                 = h2o(wth_1:wth_last),                        &
          hcn_1d                 = hcn(wth_1:wth_last),                        &
          he_1d                  = he(wth_1:wth_last),                         &
          k_1d                   = potassium(wth_1:wth_last),                  &
          li_1d                  = li(wth_1:wth_last),                         &
          n2_1d                  = n2(wth_1:wth_last),                         &
          n2o_1d                 = n2o(wth_1:wth_last),                        &
          na_1d                  = na(wth_1:wth_last),                         &
          nh3_1d                 = nh3(wth_1:wth_last),                        &
          o2_1d                  = o2(wth_1:wth_last),                         &
          o3_1d                  = o3(wth_1:wth_last),                         &
          rb_1d                  = rb(wth_1:wth_last),                         &
          so2_1d                 = so2(wth_1:wth_last),                        &
          tio_1d                 = tio(wth_1:wth_last),                        &
          vo_1d                  = vo(wth_1:wth_last),                         &
          cfc11_mix_ratio        = cfc11_mix_ratio_now,                        &
          cfc113_mix_ratio       = cfc113_mix_ratio_now,                       &
          cfc12_mix_ratio        = cfc12_mix_ratio_now,                        &
          ch4_mix_ratio          = ch4_mix_ratio_now,                          &
          co_mix_ratio           = co_mix_ratio_now,                           &
          co2_mix_ratio          = co2_mix_ratio_now,                          &
          cs_mix_ratio           = cs_mix_ratio_now,                           &
          h2_mix_ratio           = h2_mix_ratio_now,                           &
          h2o_mix_ratio          = h2o_mix_ratio_now,                          &
          hcfc22_mix_ratio       = hcfc22_mix_ratio_now,                       &
          hcn_mix_ratio          = hcn_mix_ratio_now,                          &
          he_mix_ratio           = he_mix_ratio_now,                           &
          hfc134a_mix_ratio      = hfc134a_mix_ratio_now,                      &
          k_mix_ratio            = k_mix_ratio_now,                            &
          li_mix_ratio           = li_mix_ratio_now,                           &
          n2_mix_ratio           = n2_mix_ratio_now,                           &
          n2o_mix_ratio          = n2o_mix_ratio_now,                          &
          na_mix_ratio           = na_mix_ratio_now,                           &
          nh3_mix_ratio          = nh3_mix_ratio_now,                          &
          o2_mix_ratio           = o2_mix_ratio_now,                           &
          o3_mix_ratio           = o3_mix_ratio_now,                           &
          rb_mix_ratio           = rb_mix_ratio_now,                           &
          so2_mix_ratio          = so2_mix_ratio_now,                          &
          tio_mix_ratio          = tio_mix_ratio_now,                          &
          vo_mix_ratio           = vo_mix_ratio_now,                           &
          l_ch4_well_mixed       = ch4_well_mixed,                             &
          l_co_well_mixed        = co_well_mixed,                              &
          l_co2_well_mixed       = co2_well_mixed,                             &
          l_cs_well_mixed        = cs_well_mixed,                              &
          l_h2_well_mixed        = h2_well_mixed,                              &
          l_h2o_well_mixed       = h2o_well_mixed,                             &
          l_hcn_well_mixed       = hcn_well_mixed,                             &
          l_he_well_mixed        = he_well_mixed,                              &
          l_k_well_mixed         = k_well_mixed,                               &
          l_li_well_mixed        = li_well_mixed,                              &
          l_n2_well_mixed        = n2_well_mixed,                              &
          l_n2o_well_mixed       = n2o_well_mixed,                             &
          l_na_well_mixed        = na_well_mixed,                              &
          l_nh3_well_mixed       = nh3_well_mixed,                             &
          l_o2_well_mixed        = o2_well_mixed,                              &
          l_o3_well_mixed        = o3_well_mixed,                              &
          l_rb_well_mixed        = rb_well_mixed,                              &
          l_so2_well_mixed       = so2_well_mixed,                             &
          l_tio_well_mixed       = tio_well_mixed,                             &
          l_vo_well_mixed        = vo_well_mixed,                              &
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
          layer_heat_capacity_1d = layer_heat_capacity(wth_1:wth_last),        &
          l_rayleigh             = l_rayleigh_sw,                              &
          l_mixing_ratio         = .true.,                                     &
          i_cloud_representation = i_cloud_representation,                     &
          i_overlap              = i_overlap,                                  &
          i_inhom                = i_inhom,                                    &
          i_cloud_entrapment     = i_cloud_entrapment,                         &
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
      end do
      !$OMP end do
      !$OMP end PARALLEL
    end if
  end do

end subroutine sw_code
end module sw_kernel_mod
