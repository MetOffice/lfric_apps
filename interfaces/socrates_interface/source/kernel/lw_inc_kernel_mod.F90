!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to Socrates for simple longwave increments

module lw_inc_kernel_mod

use argument_mod,      only : arg_type, &
                              GH_FIELD, GH_SCALAR, &
                              GH_REAL, GH_INTEGER, GH_LOGICAL, &
                              GH_READ, GH_READWRITE, &
                              DOMAIN, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_6, &
                              ANY_DISCONTINUOUS_SPACE_7
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: lw_inc_kernel_type
  private
  type(arg_type) :: meta_args(58) = (/ &
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    & ! lw_heating_rate_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    & ! lw_heating_rate_rtsi
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf_rtsi
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_surf_rtsi
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_toa_rtsi
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile_rtsi
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! rho_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! pressure_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! temperature_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_6), & ! t_layer_boundaries
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! d_mass
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! layer_heat_capacity
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
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_7), & ! tile_lwinc_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sulphuric
    arg_type(GH_SCALAR, GH_LOGICAL, GH_READ                                )  & ! rad_this_tstep
    /)
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: lw_inc_code
end type

public :: lw_inc_code

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                  Number of layers
!> @param[in]     n_profile                Number of columns
!> @param[in,out] lw_heating_rate_rts      LW heating rate
!> @param[in,out] lw_down_surf_rts         LW downward surface flux
!> @param[in,out] lw_up_surf_rts           LW upward surface flux
!> @param[in,out] lw_up_toa_rts            LW upward top-of-atmosphere flux
!> @param[in,out] lw_up_tile_rts           LW upward tiled surface flux
!> @param[in,out] lw_heating_rate_rtsi     LWINC heating rate
!> @param[in,out] lw_down_surf_rtsi        LWINC downward surface flux
!> @param[in,out] lw_up_surf_rtsi          LWINC upward surface flux
!> @param[in,out] lw_up_toa_rtsi           LWINC upward top-of-atmosphere flux
!> @param[in,out] lw_up_tile_rtsi          LWINC upward tiled surface flux
!> @param[in]     rho_in_wth               Density in potential temperature space
!> @param[in]     pressure_in_wth          Pressure in wth space
!> @param[in]     temperature_in_wth       Temperature in wth space
!> @param[in]     t_layer_boundaries       Temperature on radiation levels
!> @param[in]     d_mass                   Mass per square metre of radiation layers
!> @param[in]     layer_heat_capacity      Heat capacity of radiation layers
!> @param[in]     h2o                      Water vapour
!> @param[in]     co2                      Carbon dioxide
!> @param[in]     o3                       Ozone
!> @param[in]     n2o                      Dinitrogen oxide
!> @param[in]     co                       Carbon monoxide
!> @param[in]     ch4                      Methane
!> @param[in]     o2                       Oxygen
!> @param[in]     so2                      Sulphur dioxide
!> @param[in]     nh3                      Ammonia
!> @param[in]     n2                       Nitrogen
!> @param[in]     h2                       Hydrogen
!> @param[in]     he                       Helium
!> @param[in]     hcn                      Hydrogen cyanide
!> @param[in]     cs                       Cesium
!> @param[in]     potassium                Potassium
!> @param[in]     li                       Lithium
!> @param[in]     na                       Sodium
!> @param[in]     rb                       Rubidium
!> @param[in]     tio                      Titanium oxide
!> @param[in]     vo                       Vanadium oxide
!> @param[in]     mcl                      Cloud liquid field
!> @param[in]     mcf                      Cloud ice field
!> @param[in]     n_ice                    Ice number concentration
!> @param[in]     conv_frozen_number       Convective Ice number concentration
!> @param[in]     conv_liquid_mmr          Convective liquid gridbox MMR
!> @param[in]     conv_frozen_mmr          Convective frozen gridbox MMR
!> @param[in]     radiative_cloud_fraction Large scale cloud fraction
!> @param[in]     radiative_conv_fraction  Convective cloud fraction
!> @param[in]     liquid_fraction          Liquid cloud fraction field
!> @param[in]     frozen_fraction          Frozen cloud fraction field
!> @param[in]     conv_liquid_fraction     Convective liquid cloud fraction
!> @param[in]     conv_frozen_fraction     Convective frozen cloud fraction
!> @param[in]     sigma_ml                 Fractional standard deviation of liquid condensate
!> @param[in]     sigma_mi                 Fractional standard deviation of ice condensate
!> @param[in]     cloud_drop_no_conc       Cloud Droplet Number Concentration
!> @param[in]     rand_seed                Random seed field for cloud generator
!> @param[in]     n_cloud_layer            Number of cloud layers
!> @param[in]     tile_fraction            Surface tile fractions
!> @param[in]     tile_temperature         Surface tile temperature
!> @param[in]     tile_lwinc_albedo        LWINC tile albedos
!> @param[in]     sulphuric                Sulphuric acid aerosol
!> @param[in]     rad_this_tstep           Full radiation call this timestep
!> @param[in]     ndf_wth                  No. DOFs per cell for wth space
!> @param[in]     undf_wth                 No. unique of DOFs for wth space
!> @param[in]     map_wth                  Dofmap for wth space column base cell
!> @param[in]     ndf_2d                   No. of DOFs per cell for 2D space
!> @param[in]     undf_2d                  No. unique of DOFs for 2D space
!> @param[in]     map_2d                   Dofmap for 2D space column base cell
!> @param[in]     ndf_tile                 Number of DOFs per cell for tiles
!> @param[in]     undf_tile                Number of total DOFs for tiles
!> @param[in]     map_tile                 Dofmap for tile space column base cell
!> @param[in]     ndf_flux                 No. of DOFs per cell for flux space
!> @param[in]     undf_flux                No. unique of DOFs for flux space
!> @param[in]     map_flux                 Dofmap for flux space column base cell
!> @param[in]     ndf_itile                No. of DOFs per cell for itile space
!> @param[in]     undf_itile               No. unique of DOFs for itile space
!> @param[in]     map_itile                Dofmap for itile space column base cell
subroutine lw_inc_code(nlayers, n_profile,                                     &
                   lw_heating_rate_rts, lw_down_surf_rts, lw_up_surf_rts,      &
                   lw_up_toa_rts, lw_up_tile_rts,                              &
                   lw_heating_rate_rtsi, lw_down_surf_rtsi, lw_up_surf_rtsi,   &
                   lw_up_toa_rtsi, lw_up_tile_rtsi,                            &
                   rho_in_wth, pressure_in_wth, temperature_in_wth,            &
                   t_layer_boundaries, d_mass, layer_heat_capacity,            &
                   h2o, co2, o3, n2o, co, ch4, o2, so2, nh3, n2, h2, he, hcn,  &
                   cs, potassium, li, na, rb, tio, vo,                         &
                   mcl, mcf, n_ice, conv_frozen_number,                        &
                   conv_liquid_mmr, conv_frozen_mmr,                           &
                   radiative_cloud_fraction, radiative_conv_fraction,          &
                   liquid_fraction, frozen_fraction,                           &
                   conv_liquid_fraction, conv_frozen_fraction,                 &
                   sigma_ml, sigma_mi, cloud_drop_no_conc,                     &
                   rand_seed, n_cloud_layer,                                   &
                   tile_fraction, tile_temperature, tile_lwinc_albedo,         &
                   sulphuric,                                                  &
                   rad_this_tstep,                                             &
                   ndf_wth, undf_wth, map_wth,                                 &
                   ndf_2d, undf_2d, map_2d,                                    &
                   ndf_tile, undf_tile, map_tile,                              &
                   ndf_flux, undf_flux, map_flux,                              &
                   ndf_itile, undf_itile, map_itile)

  use radiation_config_mod, only: &
    i_cloud_ice_type_lwinc, i_cloud_liq_type_lwinc, &
    cloud_vertical_decorr, constant_droplet_effective_radius, &
    liu_aparam, liu_bparam
  use aerosol_config_mod, only: sulphuric_strat_climatology
  use jules_control_init_mod, only: n_surf_tile
  use socrates_init_mod, only: n_lwinc_band, &
    i_scatter_method_lwinc, &
    i_cloud_representation, i_overlap, i_inhom_inc, i_cloud_entrapment, &
    i_drop_re
  use socrates_runes, only: runes, StrDiag, ip_source_thermal
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

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, n_profile
  integer(i_def), intent(in) :: ndf_wth, undf_wth
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: ndf_flux, undf_flux
  integer(i_def), intent(in) :: ndf_itile, undf_itile

  integer(i_def), dimension(ndf_wth, n_profile),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d, n_profile),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile, n_profile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_flux, n_profile),  intent(in) :: map_flux
  integer(i_def), dimension(ndf_itile, n_profile), intent(in) :: map_itile

  real(r_def), dimension(undf_2d),   intent(inout) :: &
    lw_down_surf_rts, lw_up_surf_rts, lw_up_toa_rts
  real(r_def), dimension(undf_wth),  intent(inout) :: &
    lw_heating_rate_rts
  real(r_def), dimension(undf_tile), intent(inout) :: lw_up_tile_rts

  real(r_def), dimension(undf_2d),   intent(inout), target :: &
    lw_down_surf_rtsi, lw_up_surf_rtsi, lw_up_toa_rtsi
  real(r_def), dimension(undf_wth),  intent(inout), target :: &
    lw_heating_rate_rtsi
  real(r_def), dimension(undf_tile), intent(inout), target :: lw_up_tile_rtsi

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
  real(r_def), dimension(undf_flux), intent(in) :: t_layer_boundaries

  integer(i_def), dimension(undf_2d), intent(in) :: rand_seed, n_cloud_layer

  real(r_def), dimension(undf_tile),  intent(in) :: tile_fraction
  real(r_def), dimension(undf_tile),  intent(in) :: tile_temperature
  real(r_def), dimension(undf_itile), intent(in) :: tile_lwinc_albedo
  real(r_def), dimension(undf_wth),   intent(in) :: sulphuric

  logical(l_def), intent(in) :: rad_this_tstep

  ! Local variables for the kernel
  integer(i_def) :: n_profile_list
  integer(i_def), allocatable :: profile_list(:)
  integer(i_def) :: j, jj, k, kk, l, ll
  integer(i_def) :: wth_0, wth_1, wth_last
  integer(i_def) :: tile_1, tile_last
  integer(i_def) :: itile_1, itile_last
  integer(i_def) :: flux_0, flux_last, twod_1, twod_last
  type(StrDiag) :: lwinc_diag


  ! Set indexing
  wth_0 = map_wth(1,1)
  wth_1 = map_wth(1,1)+1
  wth_last = map_wth(1,1)+n_profile*(nlayers+1)-1
  tile_1 = map_tile(1,1)
  tile_last = map_tile(1,1)+n_profile*n_surf_tile-1
  itile_1 = map_itile(1,1)
  itile_last = map_itile(1,1)+n_profile*n_lwinc_band*n_surf_tile-1
  flux_0 = map_flux(1,1)
  flux_last = map_flux(1,1)+n_profile*(nlayers+1)-1
  twod_1 = map_2d(1,1)
  twod_last = map_2d(1,1)+n_profile-1

  ! Set pointers for diagnostic output
  if (rad_this_tstep) then
    lwinc_diag%heating_rate(0:nlayers, 1:n_profile) &
                    => lw_heating_rate_rtsi(wth_0:wth_last)
    lwinc_diag%flux_up_tile(1:n_surf_tile, 1:n_profile) &
                    => lw_up_tile_rtsi(tile_1:tile_last)
    lwinc_diag%flux_down_surf(1:n_profile) &
                    => lw_down_surf_rtsi(twod_1:twod_last)
    lwinc_diag%flux_up_surf(1:n_profile) &
                    => lw_up_surf_rtsi(twod_1:twod_last)
    lwinc_diag%flux_up_toa(1:n_profile) &
                    => lw_up_toa_rtsi(twod_1:twod_last)
  else
    allocate( lwinc_diag%heating_rate(0:nlayers, 1:n_profile) )
    allocate( lwinc_diag%flux_up_tile(1:n_surf_tile, 1:n_profile) )
    allocate( lwinc_diag%flux_down_surf(1:n_profile) )
    allocate( lwinc_diag%flux_up_surf(1:n_profile) )
    allocate( lwinc_diag%flux_up_toa(1:n_profile) )
  end if

  do k=0, nlayers
    profile_list = pack( [(l, l=1, n_profile)], &
                         n_cloud_layer(twod_1:twod_last) == k )
    n_profile_list = size(profile_list)
    if (n_profile_list > 0) then
      ! Calculate the LW increment fluxes
      call runes(n_profile_list, nlayers, lwinc_diag,                        &
        spectrum_name          = 'lwinc',                                    &
        i_source               = ip_source_thermal,                          &
        profile_list           = profile_list,                               &
        n_layer_stride         = nlayers+1,                                  &
        n_cloud_layer          = k,                                          &
        p_layer_1d             = pressure_in_wth(wth_1:wth_last),            &
        t_layer_1d             = temperature_in_wth(wth_1:wth_last),         &
        mass_1d                = d_mass(wth_1:wth_last),                     &
        density_1d             = rho_in_wth(wth_1:wth_last),                 &
        t_level_1d             = t_layer_boundaries(flux_0:flux_last),       &
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
        n_tile                 = n_surf_tile,                                &
        frac_tile_1d           = tile_fraction(tile_1:tile_last),            &
        t_tile_1d              = tile_temperature(tile_1:tile_last),         &
        albedo_diff_tile_1d    = tile_lwinc_albedo(itile_1:itile_last),      &
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
        l_mixing_ratio         = .true.,                                     &
        i_scatter_method       = i_scatter_method_lwinc,                     &
        i_cloud_representation = i_cloud_representation,                     &
        i_overlap              = i_overlap,                                  &
        i_inhom                = i_inhom_inc,                                &
        i_cloud_entrapment     = i_cloud_entrapment,                         &
        i_drop_re              = i_drop_re,                                  &
        i_st_water             = i_cloud_liq_type_lwinc,                     &
        i_st_ice               = i_cloud_ice_type_lwinc,                     &
        i_cnv_water            = i_cloud_liq_type_lwinc,                     &
        i_cnv_ice              = i_cloud_ice_type_lwinc,                     &
        l_sulphuric            = sulphuric_strat_climatology,                &
        sulphuric_1d           = sulphuric(wth_1:wth_last),                  &
        l_invert               = .true.,                                     &
        l_profile_last         = .true.)
    end if
  end do

  if (.not.rad_this_tstep) then
    ! Apply increments to radiative timestep fluxes and heating rates
    ll = twod_1 - 1
    kk = wth_0 - 1
    jj = tile_1 - 1
    do l=1, n_profile
      ll = ll + 1

      do k=0, nlayers
        kk = kk + 1
        lw_heating_rate_rts(kk) = lw_heating_rate_rts(kk) &
          - lw_heating_rate_rtsi(kk) + lwinc_diag%heating_rate(k, l)
      end do

      do j=1, n_surf_tile
        jj = jj + 1
        lw_up_tile_rts(jj) = max( 0.0_r_def, lw_up_tile_rts(jj) &
          - lw_up_tile_rtsi(jj) + lwinc_diag%flux_up_tile(j, l) )
      end do

      lw_down_surf_rts(ll) = max( 0.0_r_def, lw_down_surf_rts(ll) &
        - lw_down_surf_rtsi(ll) + lwinc_diag%flux_down_surf(l) )

      lw_up_surf_rts(ll) = max( 0.0_r_def, lw_up_surf_rts(ll) &
        - lw_up_surf_rtsi(ll) + lwinc_diag%flux_up_surf(l) )

      lw_up_toa_rts(ll) = max( 0.0_r_def, lw_up_toa_rts(ll) &
        - lw_up_toa_rtsi(ll) + lwinc_diag%flux_up_toa(l) )
    end do
    deallocate( lwinc_diag%flux_up_toa )
    deallocate( lwinc_diag%flux_up_surf )
    deallocate( lwinc_diag%flux_down_surf )
    deallocate( lwinc_diag%flux_up_tile )
    deallocate( lwinc_diag%heating_rate )
  end if

end subroutine lw_inc_code
end module lw_inc_kernel_mod
