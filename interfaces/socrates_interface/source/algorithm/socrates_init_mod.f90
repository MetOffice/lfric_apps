!------------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
! @brief Initialisation for the Socrates radiation code

module socrates_init_mod

use constants_mod, only: r_def, i_def, l_def

implicit none

integer(i_def), pointer :: n_band_exclude(:)
integer(i_def), pointer :: index_exclude(:,:)
real(r_def), pointer :: wavelength_short(:)
real(r_def), pointer :: wavelength_long(:)
real(r_def), pointer :: weight_blue(:)

integer(i_def) :: n_sw_band, n_lw_band
integer(i_def), allocatable, target :: sw_n_band_exclude(:)
integer(i_def), allocatable, target :: sw_index_exclude(:, :)
integer(i_def), allocatable, target :: lw_n_band_exclude(:)
integer(i_def), allocatable, target :: lw_index_exclude(:, :)
real(r_def), allocatable, target :: &
  sw_wavelength_short(:), sw_wavelength_long(:), sw_weight_blue(:)
real(r_def), allocatable, target :: &
  lw_wavelength_short(:), lw_wavelength_long(:)

integer(i_def) :: n_swinc_band, n_lwinc_band
integer(i_def), allocatable, target :: swinc_n_band_exclude(:)
integer(i_def), allocatable, target :: swinc_index_exclude(:, :)
integer(i_def), allocatable, target :: lwinc_n_band_exclude(:)
integer(i_def), allocatable, target :: lwinc_index_exclude(:, :)
real(r_def), allocatable, target :: &
  swinc_wavelength_short(:), swinc_wavelength_long(:), swinc_weight_blue(:)
real(r_def), allocatable, target :: &
  lwinc_wavelength_short(:), lwinc_wavelength_long(:)

integer(i_def) :: n_pathway, sw_550nm, lw_11um

integer(i_def) :: i_scatter_method_lw, i_scatter_method_lwinc
integer(i_def) :: i_cloud_representation, i_overlap, i_inhom, i_inhom_inc
integer(i_def) :: i_cloud_entrapment
integer(i_def) :: i_drop_re
logical(l_def) :: l_orog

private
public :: socrates_init, &
  wavelength_short, wavelength_long, &
  weight_blue, n_band_exclude, index_exclude, &
  n_sw_band, sw_n_band_exclude, sw_index_exclude, &
  sw_wavelength_short, sw_wavelength_long, sw_weight_blue, &
  n_lw_band, lw_n_band_exclude, lw_index_exclude, &
  lw_wavelength_short, lw_wavelength_long, &
  n_swinc_band, swinc_n_band_exclude, swinc_index_exclude, &
  swinc_wavelength_short, swinc_wavelength_long, swinc_weight_blue, &
  n_lwinc_band, lwinc_n_band_exclude, lwinc_index_exclude, &
  lwinc_wavelength_short, lwinc_wavelength_long, &
  n_pathway, sw_550nm, lw_11um, &
  i_scatter_method_lw, i_scatter_method_lwinc, &
  i_cloud_representation, i_overlap, i_inhom, i_inhom_inc, &
  i_cloud_entrapment, i_drop_re, l_orog

contains

subroutine socrates_init()

  use cosp_config_mod, only: l_cosp
  use radiation_config_mod, only: &
    spectral_file_sw, spectral_file_lw, mcica_data_file, &
    scatter_method_lw, &
    scatter_method_lw_full, scatter_method_lw_none, &
    scatter_method_lw_approx, scatter_method_lw_hybrid, &
    scatter_method_lwinc, &
    scatter_method_lwinc_full, scatter_method_lwinc_none, &
    scatter_method_lwinc_approx, scatter_method_lwinc_hybrid, &
    cloud_representation, &
    cloud_representation_no_cloud, &
    cloud_representation_liquid_and_ice, &
    cloud_representation_combined, &
    cloud_representation_conv_strat_liq_ice, &
    cloud_representation_split, &
    cloud_overlap, &
    cloud_overlap_maximum_random, &
    cloud_overlap_random, &
    cloud_overlap_exponential_random, &
    cloud_inhomogeneity, &
    cloud_inhomogeneity_homogeneous, &
    cloud_inhomogeneity_scaling, &
    cloud_inhomogeneity_mcica, &
    cloud_inhomogeneity_cairns, &
    cloud_inhomogeneity_tripleclouds, &
    cloud_entrapment, &
    cloud_entrapment_max, &
    cloud_entrapment_zero, &
    droplet_effective_radius, &
    droplet_effective_radius_constant, &
    droplet_effective_radius_liu, &
    l_inc_radstep, spectral_file_swinc, spectral_file_lwinc, &
    topography, topography_flat
  use radiative_gases_config_mod, only: &
    cfc11_rad_opt, cfc11_rad_opt_off, &
    cfc113_rad_opt, cfc113_rad_opt_off, &
    cfc12_rad_opt, cfc12_rad_opt_off, &
    ch4_rad_opt, ch4_rad_opt_off, &
    co_rad_opt, co_rad_opt_off, &
    co2_rad_opt, co2_rad_opt_off, &
    cs_rad_opt, cs_rad_opt_off, &
    h2_rad_opt, h2_rad_opt_off, &
    h2o_rad_opt, h2o_rad_opt_off, &
    hcfc22_rad_opt, hcfc22_rad_opt_off, &
    hcn_rad_opt, hcn_rad_opt_off, &
    he_rad_opt, he_rad_opt_off, &
    hfc134a_rad_opt, hfc134a_rad_opt_off, &
    k_rad_opt, k_rad_opt_off, &
    li_rad_opt, li_rad_opt_off, &
    n2_rad_opt, n2_rad_opt_off, &
    n2o_rad_opt, n2o_rad_opt_off, &
    na_rad_opt, na_rad_opt_off, &
    nh3_rad_opt, nh3_rad_opt_off, &
    o2_rad_opt, o2_rad_opt_off, &
    o3_rad_opt, o3_rad_opt_off, &
    rb_rad_opt, rb_rad_opt_off, &
    so2_rad_opt, so2_rad_opt_off, &
    tio_rad_opt, tio_rad_opt_off, &
    vo_rad_opt, vo_rad_opt_off
  use rad_ccf, only: set_socrates_constants
  use socrates_runes, only: &
    ip_scatter_full, ip_scatter_none, ip_scatter_approx, ip_scatter_hybrid, &
    ip_cloud_representation_off, ip_cloud_representation_ice_water, &
    ip_cloud_representation_combine_ice_water, &
    ip_cloud_representation_csiw, ip_cloud_representation_split_ice_water, &
    ip_overlap_max_random, ip_overlap_random, ip_overlap_exponential_random, &
    ip_inhom_homogeneous, ip_inhom_scaling, ip_inhom_mcica, ip_inhom_cairns, &
    ip_inhom_tripleclouds_2019, &
    ip_cloud_entrapment_max, ip_cloud_entrapment_zero, &
    ip_droplet_re_default, ip_droplet_re_constant, ip_droplet_re_liu
  use socrates_set_spectrum, only: &
    set_spectrum, get_spectrum, get_spectral_band, set_mcica
  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_INFO

  implicit none

  integer(i_def) :: i_band
  integer(i_def) :: n_spectral_files

  ! Switches for radiatively active gases
  logical(l_def) :: l_cfc11
  logical(l_def) :: l_cfc113
  logical(l_def) :: l_cfc12
  logical(l_def) :: l_ch4
  logical(l_def) :: l_co
  logical(l_def) :: l_co2
  logical(l_def) :: l_cs
  logical(l_def) :: l_h2
  logical(l_def) :: l_h2o
  logical(l_def) :: l_hcfc22
  logical(l_def) :: l_hcn
  logical(l_def) :: l_he
  logical(l_def) :: l_hfc134a
  logical(l_def) :: l_k
  logical(l_def) :: l_li
  logical(l_def) :: l_n2
  logical(l_def) :: l_n2o
  logical(l_def) :: l_na
  logical(l_def) :: l_nh3
  logical(l_def) :: l_o2
  logical(l_def) :: l_o3
  logical(l_def) :: l_rb
  logical(l_def) :: l_so2
  logical(l_def) :: l_tio
  logical(l_def) :: l_vo
  l_cfc11 = (cfc11_rad_opt /= cfc11_rad_opt_off)
  l_cfc113 = (cfc113_rad_opt /= cfc113_rad_opt_off)
  l_cfc12 = (cfc12_rad_opt /= cfc12_rad_opt_off)
  l_ch4 = (ch4_rad_opt /= ch4_rad_opt_off)
  l_co = (co_rad_opt /= co_rad_opt_off)
  l_co2 = (co2_rad_opt /= co2_rad_opt_off)
  l_cs = (cs_rad_opt /= cs_rad_opt_off)
  l_h2 = (h2_rad_opt /= h2_rad_opt_off)
  l_h2o = (h2o_rad_opt /= h2o_rad_opt_off)
  l_hcfc22 = (hcfc22_rad_opt /= hcfc22_rad_opt_off)
  l_hcn = (hcn_rad_opt /= hcn_rad_opt_off)
  l_he = (he_rad_opt /= he_rad_opt_off)
  l_hfc134a = (hfc134a_rad_opt /= hfc134a_rad_opt_off)
  l_k = (k_rad_opt /= k_rad_opt_off)
  l_li = (li_rad_opt /= li_rad_opt_off)
  l_n2 = (n2_rad_opt /= n2_rad_opt_off)
  l_n2o = (n2o_rad_opt /= n2o_rad_opt_off)
  l_na = (na_rad_opt /= na_rad_opt_off)
  l_nh3 = (nh3_rad_opt /= nh3_rad_opt_off)
  l_o2 = (o2_rad_opt /= o2_rad_opt_off)
  l_o3 = (o3_rad_opt /= o3_rad_opt_off)
  l_rb = (rb_rad_opt /= rb_rad_opt_off)
  l_so2 = (so2_rad_opt /= so2_rad_opt_off)
  l_tio = (tio_rad_opt /= tio_rad_opt_off)
  l_vo = (vo_rad_opt /= vo_rad_opt_off)

  ! Set constants in the socrates modules
  call set_socrates_constants()

  if (l_inc_radstep) then
    ! Read spectral files for incremental radiative timestepping
    n_spectral_files = 4
    call set_spectrum(                           &
      n_instances      = n_spectral_files,       &
      spectrum_name    = 'swinc',                &
      spectral_file    = spectral_file_swinc,    &
      l_cfc11          = l_cfc11,                &
      l_cfc113         = l_cfc113,               &
      l_cfc12          = l_cfc12,                &
      l_ch4            = l_ch4,                  &
      l_co             = l_co,                   &
      l_co2            = l_co2,                  &
      l_cs             = l_cs,                   &
      l_h2             = l_h2,                   &
      l_h2o            = l_h2o,                  &
      l_hcfc22         = l_hcfc22,               &
      l_hcn            = l_hcn,                  &
      l_he             = l_he,                   &
      l_hfc134a        = l_hfc134a,              &
      l_k              = l_k,                    &
      l_li             = l_li,                   &
      l_n2             = l_n2,                   &
      l_n2o            = l_n2o,                  &
      l_na             = l_na,                   &
      l_nh3            = l_nh3,                  &
      l_o2             = l_o2,                   &
      l_o3             = l_o3,                   &
      l_rb             = l_rb,                   &
      l_so2            = l_so2,                  &
      l_tio            = l_tio,                  &
      l_vo             = l_vo                    )
    call get_spectrum(                           &
      spectrum_name    = 'swinc',                &
      n_band           = n_swinc_band,           &
      n_band_exclude   = swinc_n_band_exclude,   &
      index_exclude    = swinc_index_exclude,    &
      wavelength_short = swinc_wavelength_short, &
      wavelength_long  = swinc_wavelength_long,  &
      weight_blue      = swinc_weight_blue )
    call set_spectrum(                           &
      spectrum_name    = 'lwinc',                &
      spectral_file    = spectral_file_lwinc,    &
      l_cfc11          = l_cfc11,                &
      l_cfc113         = l_cfc113,               &
      l_cfc12          = l_cfc12,                &
      l_ch4            = l_ch4,                  &
      l_co             = l_co,                   &
      l_co2            = l_co2,                  &
      l_cs             = l_cs,                   &
      l_h2             = l_h2,                   &
      l_h2o            = l_h2o,                  &
      l_hcfc22         = l_hcfc22,               &
      l_hcn            = l_hcn,                  &
      l_he             = l_he,                   &
      l_hfc134a        = l_hfc134a,              &
      l_k              = l_k,                    &
      l_li             = l_li,                   &
      l_n2             = l_n2,                   &
      l_n2o            = l_n2o,                  &
      l_na             = l_na,                   &
      l_nh3            = l_nh3,                  &
      l_o2             = l_o2,                   &
      l_o3             = l_o3,                   &
      l_rb             = l_rb,                   &
      l_so2            = l_so2,                  &
      l_tio            = l_tio,                  &
      l_vo             = l_vo                    )
    call get_spectrum(                           &
      spectrum_name    = 'lwinc',                &
      n_band           = n_lwinc_band,           &
      n_band_exclude   = lwinc_n_band_exclude,   &
      index_exclude    = lwinc_index_exclude,    &
      wavelength_short = lwinc_wavelength_short, &
      wavelength_long  = lwinc_wavelength_long )
  end if

  ! Read spectral files
  call set_spectrum(                        &
    spectrum_name    = 'sw',                &
    spectral_file    = spectral_file_sw,    &
    l_cfc11          = l_cfc11,             &
    l_cfc113         = l_cfc113,            &
    l_cfc12          = l_cfc12,             &
    l_ch4            = l_ch4,               &
    l_co             = l_co,                &
    l_co2            = l_co2,               &
    l_cs             = l_cs,                &
    l_h2             = l_h2,                &
    l_h2o            = l_h2o,               &
    l_hcfc22         = l_hcfc22,            &
    l_hcn            = l_hcn,               &
    l_he             = l_he,                &
    l_hfc134a        = l_hfc134a,           &
    l_k              = l_k,                 &
    l_li             = l_li,                &
    l_n2             = l_n2,                &
    l_n2o            = l_n2o,               &
    l_na             = l_na,                &
    l_nh3            = l_nh3,               &
    l_o2             = l_o2,                &
    l_o3             = l_o3,                &
    l_rb             = l_rb,                &
    l_so2            = l_so2,               &
    l_tio            = l_tio,               &
    l_vo             = l_vo                 )
  call get_spectrum(                        &
    spectrum_name    = 'sw',                &
    n_band           = n_sw_band,           &
    n_band_exclude   = sw_n_band_exclude,   &
    index_exclude    = sw_index_exclude,    &
    wavelength_short = sw_wavelength_short, &
    wavelength_long  = sw_wavelength_long,  &
    weight_blue      = sw_weight_blue,      &
    n_pathway        = n_pathway )

  sw_550nm = get_spectral_band(             &
    wavelength       = 550.0e-9_r_def,      &
    spectrum_name    = 'sw' )

  call log_event( 'SW bands:', LOG_LEVEL_INFO )
  do i_band=1, n_sw_band
    write( log_scratch_space, '(I3,2E16.8)' ) &
      i_band, sw_wavelength_short(i_band), sw_wavelength_long(i_band)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

  call set_spectrum(                        &
    spectrum_name    = 'lw',                &
    spectral_file    = spectral_file_lw,    &
    l_cfc11          = l_cfc11,             &
    l_cfc113         = l_cfc113,            &
    l_cfc12          = l_cfc12,             &
    l_ch4            = l_ch4,               &
    l_co             = l_co,                &
    l_co2            = l_co2,               &
    l_cs             = l_cs,                &
    l_h2             = l_h2,                &
    l_h2o            = l_h2o,               &
    l_hcfc22         = l_hcfc22,            &
    l_hcn            = l_hcn,               &
    l_he             = l_he,                &
    l_hfc134a        = l_hfc134a,           &
    l_k              = l_k,                 &
    l_li             = l_li,                &
    l_n2             = l_n2,                &
    l_n2o            = l_n2o,               &
    l_na             = l_na,                &
    l_nh3            = l_nh3,               &
    l_o2             = l_o2,                &
    l_o3             = l_o3,                &
    l_rb             = l_rb,                &
    l_so2            = l_so2,               &
    l_tio            = l_tio,               &
    l_vo             = l_vo                 )
  call get_spectrum(                        &
    spectrum_name    = 'lw',                &
    n_band           = n_lw_band,           &
    n_band_exclude   = lw_n_band_exclude,   &
    index_exclude    = lw_index_exclude,    &
    wavelength_short = lw_wavelength_short, &
    wavelength_long  = lw_wavelength_long )

  lw_11um = get_spectral_band(              &
    wavelength       = 11.0e-6_r_def,       &
    spectrum_name    = 'lw' )

  call log_event( 'LW bands:', LOG_LEVEL_INFO )
  do i_band=1, n_lw_band
    write( log_scratch_space, '(I3,2E16.8)' ) &
      i_band, lw_wavelength_short(i_band), lw_wavelength_long(i_band)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

  if ( (cloud_representation /= cloud_representation_no_cloud) .and. &
       ((cloud_inhomogeneity == cloud_inhomogeneity_mcica) .or. l_cosp) ) then
    ! Read MCICA data file
    call set_mcica(mcica_data_file, 'sw', 'lw')
    call log_event( 'Read MCICA data file.', LOG_LEVEL_INFO )
  end if

  ! LW scattering method
  select case (scatter_method_lw)
  case (scatter_method_lw_full)
    i_scatter_method_lw = ip_scatter_full
  case (scatter_method_lw_none)
    i_scatter_method_lw = ip_scatter_none
  case (scatter_method_lw_approx)
    i_scatter_method_lw = ip_scatter_approx
  case (scatter_method_lw_hybrid)
    i_scatter_method_lw = ip_scatter_hybrid
  case default
    i_scatter_method_lw = ip_scatter_full
  end select

  ! LW increment scattering method
  select case (scatter_method_lwinc)
  case (scatter_method_lwinc_full)
    i_scatter_method_lwinc = ip_scatter_full
  case (scatter_method_lwinc_none)
    i_scatter_method_lwinc = ip_scatter_none
  case (scatter_method_lwinc_approx)
    i_scatter_method_lwinc = ip_scatter_approx
  case (scatter_method_lwinc_hybrid)
    i_scatter_method_lwinc = ip_scatter_hybrid
  case default
    i_scatter_method_lwinc = ip_scatter_full
  end select

  ! Properties of clouds
  select case (cloud_representation)
  case (cloud_representation_no_cloud)
    i_cloud_representation = ip_cloud_representation_off
  case (cloud_representation_liquid_and_ice)
    i_cloud_representation = ip_cloud_representation_ice_water
  case (cloud_representation_combined)
    i_cloud_representation = ip_cloud_representation_combine_ice_water
  case (cloud_representation_conv_strat_liq_ice)
    i_cloud_representation = ip_cloud_representation_csiw
  case (cloud_representation_split)
    i_cloud_representation = ip_cloud_representation_split_ice_water
  case default
    i_cloud_representation = ip_cloud_representation_off
  end select

  select case (cloud_overlap)
  case (cloud_overlap_maximum_random)
    i_overlap = ip_overlap_max_random
  case (cloud_overlap_random)
    i_overlap = ip_overlap_random
  case (cloud_overlap_exponential_random)
    i_overlap = ip_overlap_exponential_random
  case default
    i_overlap = ip_overlap_max_random
  end select

  select case (cloud_inhomogeneity)
  case (cloud_inhomogeneity_homogeneous)
    i_inhom = ip_inhom_homogeneous
  case (cloud_inhomogeneity_scaling)
    i_inhom = ip_inhom_scaling
  case (cloud_inhomogeneity_mcica)
    i_inhom = ip_inhom_mcica
  case (cloud_inhomogeneity_cairns)
    i_inhom = ip_inhom_cairns
  case (cloud_inhomogeneity_tripleclouds)
    i_inhom = ip_inhom_tripleclouds_2019
  case default
    i_inhom = ip_inhom_homogeneous
  end select

  ! Increment radiation calls should not use MCICA as there will not
  ! be enough k-terms to sample the cloud field.
  if (i_inhom == ip_inhom_mcica) then
    i_inhom_inc = ip_inhom_scaling
  else
    i_inhom_inc = i_inhom
  end if

  select case (cloud_entrapment)
  case (cloud_entrapment_max)
    i_cloud_entrapment = ip_cloud_entrapment_max
  case (cloud_entrapment_zero)
    i_cloud_entrapment = ip_cloud_entrapment_zero
  case default
    i_cloud_entrapment = ip_cloud_entrapment_zero
  end select

  select case (droplet_effective_radius)
  case (droplet_effective_radius_constant)
    i_drop_re = ip_droplet_re_constant
  case (droplet_effective_radius_liu)
    i_drop_re = ip_droplet_re_liu
  case default
    i_drop_re = ip_droplet_re_default
  end select

  ! Set logical for the orographic correction
  if (topography == topography_flat) then
    l_orog = .false.
  else
    l_orog = .true.
  end if

end subroutine socrates_init
end module socrates_init_mod
