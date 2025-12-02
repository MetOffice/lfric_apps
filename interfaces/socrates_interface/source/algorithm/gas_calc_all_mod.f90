module gas_calc_all_mod

  use constants_mod, only: i_def, l_def, r_def, r_um
  use gas_calc_mod, only: gas_calc
  use radiative_gases_config_mod, only : &
    cfc11_clim_fcg_levls, cfc11_clim_fcg_nyears, &
    cfc11_clim_fcg_years, cfc11_clim_fcg_rates, &
    cfc11_mix_ratio, cfc11_rad_opt, cfc11_rad_opt_off, &
    cfc11_rad_opt_constant, cfc11_rad_opt_time_varying, &
    cfc113_clim_fcg_levls, cfc113_clim_fcg_nyears, &
    cfc113_clim_fcg_years, cfc113_clim_fcg_rates, &
    cfc113_mix_ratio, cfc113_rad_opt, cfc113_rad_opt_off, &
    cfc113_rad_opt_constant, cfc113_rad_opt_time_varying, &
    cfc12_clim_fcg_levls, cfc12_clim_fcg_nyears, &
    cfc12_clim_fcg_years, cfc12_clim_fcg_rates, &
    cfc12_mix_ratio, cfc12_rad_opt, cfc12_rad_opt_off, &
    cfc12_rad_opt_constant, cfc12_rad_opt_time_varying, &
    ch4_clim_fcg_levls, ch4_clim_fcg_nyears, &
    ch4_clim_fcg_years, ch4_clim_fcg_rates, &
    ch4_mix_ratio, ch4_rad_opt, ch4_rad_opt_off, &
    ch4_rad_opt_constant, ch4_rad_opt_time_varying, &
    ch4_rad_opt_prognostic, ch4_rad_opt_ancil, &
    co_clim_fcg_levls, co_clim_fcg_nyears, &
    co_clim_fcg_years, co_clim_fcg_rates, &
    co_mix_ratio, co_rad_opt, co_rad_opt_off, &
    co_rad_opt_constant, co_rad_opt_time_varying, &
    co_rad_opt_prognostic, co_rad_opt_ancil, &
    co2_clim_fcg_levls, co2_clim_fcg_nyears, &
    co2_clim_fcg_years, co2_clim_fcg_rates, &
    co2_mix_ratio, co2_rad_opt, co2_rad_opt_off, &
    co2_rad_opt_constant, co2_rad_opt_time_varying, &
    co2_rad_opt_prognostic, co2_rad_opt_ancil, &
    cs_clim_fcg_levls, cs_clim_fcg_nyears, &
    cs_clim_fcg_years, cs_clim_fcg_rates, &
    cs_mix_ratio, cs_rad_opt, cs_rad_opt_off, &
    cs_rad_opt_constant, cs_rad_opt_time_varying, &
    cs_rad_opt_prognostic, cs_rad_opt_ancil, &
    h2_clim_fcg_levls, h2_clim_fcg_nyears, &
    h2_clim_fcg_years, h2_clim_fcg_rates, &
    h2_mix_ratio, h2_rad_opt, h2_rad_opt_off, &
    h2_rad_opt_constant, h2_rad_opt_time_varying, &
    h2_rad_opt_prognostic, h2_rad_opt_ancil, &
    h2o_clim_fcg_levls, h2o_clim_fcg_nyears, &
    h2o_clim_fcg_years, h2o_clim_fcg_rates, &
    h2o_mix_ratio, h2o_rad_opt, h2o_rad_opt_off, &
    h2o_rad_opt_constant, h2o_rad_opt_time_varying, &
    h2o_rad_opt_prognostic, h2o_rad_opt_ancil, &
    hcfc22_clim_fcg_levls, hcfc22_clim_fcg_nyears, &
    hcfc22_clim_fcg_years, hcfc22_clim_fcg_rates, &
    hcfc22_mix_ratio, hcfc22_rad_opt, hcfc22_rad_opt_off, &
    hcfc22_rad_opt_constant, hcfc22_rad_opt_time_varying, &
    hcn_clim_fcg_levls, hcn_clim_fcg_nyears, &
    hcn_clim_fcg_years, hcn_clim_fcg_rates, &
    hcn_mix_ratio, hcn_rad_opt, hcn_rad_opt_off, &
    hcn_rad_opt_constant, hcn_rad_opt_time_varying, &
    hcn_rad_opt_prognostic, hcn_rad_opt_ancil, &
    he_clim_fcg_levls, he_clim_fcg_nyears, &
    he_clim_fcg_years, he_clim_fcg_rates, &
    he_mix_ratio, he_rad_opt, he_rad_opt_off, &
    he_rad_opt_constant, he_rad_opt_time_varying, &
    he_rad_opt_prognostic, he_rad_opt_ancil, &
    hfc134a_clim_fcg_levls, hfc134a_clim_fcg_nyears, &
    hfc134a_clim_fcg_years, hfc134a_clim_fcg_rates, &
    hfc134a_mix_ratio, hfc134a_rad_opt, hfc134a_rad_opt_off, &
    hfc134a_rad_opt_constant, hfc134a_rad_opt_time_varying, &
    k_clim_fcg_levls, k_clim_fcg_nyears, &
    k_clim_fcg_years, k_clim_fcg_rates, &
    k_mix_ratio, k_rad_opt, k_rad_opt_off, &
    k_rad_opt_constant, k_rad_opt_time_varying, &
    k_rad_opt_prognostic, k_rad_opt_ancil, &
    li_clim_fcg_levls, li_clim_fcg_nyears, &
    li_clim_fcg_years, li_clim_fcg_rates, &
    li_mix_ratio, li_rad_opt, li_rad_opt_off, &
    li_rad_opt_constant, li_rad_opt_time_varying, &
    li_rad_opt_prognostic, li_rad_opt_ancil, &
    n2_clim_fcg_levls, n2_clim_fcg_nyears, &
    n2_clim_fcg_years, n2_clim_fcg_rates, &
    n2_mix_ratio, n2_rad_opt, n2_rad_opt_off, &
    n2_rad_opt_constant, n2_rad_opt_time_varying, &
    n2_rad_opt_prognostic, n2_rad_opt_ancil, &
    n2o_clim_fcg_levls, n2o_clim_fcg_nyears, &
    n2o_clim_fcg_years, n2o_clim_fcg_rates, &
    n2o_mix_ratio, n2o_rad_opt, n2o_rad_opt_off, &
    n2o_rad_opt_constant, n2o_rad_opt_time_varying, &
    n2o_rad_opt_prognostic, n2o_rad_opt_ancil, &
    na_clim_fcg_levls, na_clim_fcg_nyears, &
    na_clim_fcg_years, na_clim_fcg_rates, &
    na_mix_ratio, na_rad_opt, na_rad_opt_off, &
    na_rad_opt_constant, na_rad_opt_time_varying, &
    na_rad_opt_prognostic, na_rad_opt_ancil, &
    nh3_clim_fcg_levls, nh3_clim_fcg_nyears, &
    nh3_clim_fcg_years, nh3_clim_fcg_rates, &
    nh3_mix_ratio, nh3_rad_opt, nh3_rad_opt_off, &
    nh3_rad_opt_constant, nh3_rad_opt_time_varying, &
    nh3_rad_opt_prognostic, nh3_rad_opt_ancil, &
    o2_clim_fcg_levls, o2_clim_fcg_nyears, &
    o2_clim_fcg_years, o2_clim_fcg_rates, &
    o2_mix_ratio, o2_rad_opt, o2_rad_opt_off, &
    o2_rad_opt_constant, o2_rad_opt_time_varying, &
    o2_rad_opt_prognostic, o2_rad_opt_ancil, &
    o3_clim_fcg_levls, o3_clim_fcg_nyears, &
    o3_clim_fcg_years, o3_clim_fcg_rates, &
    o3_mix_ratio, o3_rad_opt, o3_rad_opt_off, &
    o3_rad_opt_constant, o3_rad_opt_time_varying, &
    o3_rad_opt_profile, o3_rad_opt_prognostic, o3_rad_opt_ancil, &
    rb_clim_fcg_levls, rb_clim_fcg_nyears, &
    rb_clim_fcg_years, rb_clim_fcg_rates, &
    rb_mix_ratio, rb_rad_opt, rb_rad_opt_off, &
    rb_rad_opt_constant, rb_rad_opt_time_varying, &
    rb_rad_opt_prognostic, rb_rad_opt_ancil, &
    so2_clim_fcg_levls, so2_clim_fcg_nyears, &
    so2_clim_fcg_years, so2_clim_fcg_rates, &
    so2_mix_ratio, so2_rad_opt, so2_rad_opt_off, &
    so2_rad_opt_constant, so2_rad_opt_time_varying, &
    so2_rad_opt_prognostic, so2_rad_opt_ancil, &
    tio_clim_fcg_levls, tio_clim_fcg_nyears, &
    tio_clim_fcg_years, tio_clim_fcg_rates, &
    tio_mix_ratio, tio_rad_opt, tio_rad_opt_off, &
    tio_rad_opt_constant, tio_rad_opt_time_varying, &
    tio_rad_opt_prognostic, tio_rad_opt_ancil, &
    vo_clim_fcg_levls, vo_clim_fcg_nyears, &
    vo_clim_fcg_years, vo_clim_fcg_rates, &
    vo_mix_ratio, vo_rad_opt, vo_rad_opt_off, &
    vo_rad_opt_constant, vo_rad_opt_time_varying, &
    vo_rad_opt_prognostic, vo_rad_opt_ancil

  implicit none

  ! Flags for treating gas mixing ratios as well-mixed
  logical(l_def), public :: ch4_well_mixed
  logical(l_def), public :: co_well_mixed
  logical(l_def), public :: co2_well_mixed
  logical(l_def), public :: cs_well_mixed
  logical(l_def), public :: h2_well_mixed
  logical(l_def), public :: h2o_well_mixed
  logical(l_def), public :: hcn_well_mixed
  logical(l_def), public :: he_well_mixed
  logical(l_def), public :: k_well_mixed
  logical(l_def), public :: li_well_mixed
  logical(l_def), public :: n2_well_mixed
  logical(l_def), public :: n2o_well_mixed
  logical(l_def), public :: na_well_mixed
  logical(l_def), public :: nh3_well_mixed
  logical(l_def), public :: o2_well_mixed
  logical(l_def), public :: o3_well_mixed
  logical(l_def), public :: rb_well_mixed
  logical(l_def), public :: so2_well_mixed
  logical(l_def), public :: tio_well_mixed
  logical(l_def), public :: vo_well_mixed

  ! Well-mixed gas mixing ratios
  real(r_def), public :: cfc11_mix_ratio_now
  real(r_def), public :: cfc113_mix_ratio_now
  real(r_def), public :: cfc12_mix_ratio_now
  real(r_def), public :: ch4_mix_ratio_now
  real(r_def), public :: co_mix_ratio_now
  real(r_def), public :: co2_mix_ratio_now
  real(r_def), public :: cs_mix_ratio_now
  real(r_def), public :: h2_mix_ratio_now
  real(r_def), public :: h2o_mix_ratio_now
  real(r_def), public :: hcfc22_mix_ratio_now
  real(r_def), public :: hcn_mix_ratio_now
  real(r_def), public :: he_mix_ratio_now
  real(r_def), public :: hfc134a_mix_ratio_now
  real(r_def), public :: k_mix_ratio_now
  real(r_def), public :: li_mix_ratio_now
  real(r_def), public :: n2_mix_ratio_now
  real(r_def), public :: n2o_mix_ratio_now
  real(r_def), public :: na_mix_ratio_now
  real(r_def), public :: nh3_mix_ratio_now
  real(r_def), public :: o2_mix_ratio_now
  real(r_def), public :: o3_mix_ratio_now
  real(r_def), public :: rb_mix_ratio_now
  real(r_def), public :: so2_mix_ratio_now
  real(r_def), public :: tio_mix_ratio_now
  real(r_def), public :: vo_mix_ratio_now

contains

  subroutine gas_calc_all()

    implicit none

    ! Set MMR for CFC11
    select case ( cfc11_rad_opt )
    case ( cfc11_rad_opt_off )
      cfc11_mix_ratio_now = 0.0_r_def
    case ( cfc11_rad_opt_constant )
      cfc11_mix_ratio_now = cfc11_mix_ratio
    case ( cfc11_rad_opt_time_varying )
      call gas_calc( cfc11_mix_ratio_now,    &
                     cfc11_clim_fcg_nyears,  &
                     cfc11_clim_fcg_years,   &
                     cfc11_clim_fcg_levls,   &
                     cfc11_clim_fcg_rates )
    end select

    ! Set MMR for CFC113
    select case ( cfc113_rad_opt )
    case ( cfc113_rad_opt_off )
      cfc113_mix_ratio_now = 0.0_r_def
    case ( cfc113_rad_opt_constant )
      cfc113_mix_ratio_now = cfc113_mix_ratio
    case ( cfc113_rad_opt_time_varying )
      call gas_calc( cfc113_mix_ratio_now,    &
                     cfc113_clim_fcg_nyears,  &
                     cfc113_clim_fcg_years,   &
                     cfc113_clim_fcg_levls,   &
                     cfc113_clim_fcg_rates )
    end select

    ! Set MMR for CFC12
    select case ( cfc12_rad_opt )
    case ( cfc12_rad_opt_off )
      cfc12_mix_ratio_now = 0.0_r_def
    case ( cfc12_rad_opt_constant )
      cfc12_mix_ratio_now = cfc12_mix_ratio
    case ( cfc12_rad_opt_time_varying )
      call gas_calc( cfc12_mix_ratio_now,    &
                     cfc12_clim_fcg_nyears,  &
                     cfc12_clim_fcg_years,   &
                     cfc12_clim_fcg_levls,   &
                     cfc12_clim_fcg_rates )
    end select

    ! Set MMR for CH4
    select case ( ch4_rad_opt )
    case ( ch4_rad_opt_off )
      ch4_mix_ratio_now = 0.0_r_def
      ch4_well_mixed = .true.
    case ( ch4_rad_opt_constant )
      ch4_mix_ratio_now = ch4_mix_ratio
      ch4_well_mixed = .true.
    case ( ch4_rad_opt_time_varying )
      call gas_calc( ch4_mix_ratio_now,    &
                     ch4_clim_fcg_nyears,  &
                     ch4_clim_fcg_years,   &
                     ch4_clim_fcg_levls,   &
                     ch4_clim_fcg_rates )
      ch4_well_mixed = .true.
    case ( ch4_rad_opt_prognostic, ch4_rad_opt_ancil )
      ch4_well_mixed = .false.
    end select

    ! Set MMR for CO
    select case ( co_rad_opt )
    case ( co_rad_opt_off )
      co_mix_ratio_now = 0.0_r_def
      co_well_mixed = .true.
    case ( co_rad_opt_constant )
      co_mix_ratio_now = co_mix_ratio
      co_well_mixed = .true.
    case ( co_rad_opt_time_varying )
      call gas_calc( co_mix_ratio_now,    &
                     co_clim_fcg_nyears,  &
                     co_clim_fcg_years,   &
                     co_clim_fcg_levls,   &
                     co_clim_fcg_rates )
      co_well_mixed = .true.
    case ( co_rad_opt_prognostic, co_rad_opt_ancil )
      co_well_mixed = .false.
    end select

    ! Set MMR for CO2
    select case ( co2_rad_opt )
    case ( co2_rad_opt_off )
      co2_mix_ratio_now = 0.0_r_def
      co2_well_mixed = .true.
    case ( co2_rad_opt_constant )
      co2_mix_ratio_now = co2_mix_ratio
      co2_well_mixed = .true.
    case ( co2_rad_opt_time_varying )
      call gas_calc( co2_mix_ratio_now,    &
                     co2_clim_fcg_nyears,  &
                     co2_clim_fcg_years,   &
                     co2_clim_fcg_levls,   &
                     co2_clim_fcg_rates )
      co2_well_mixed = .true.
    case ( co2_rad_opt_prognostic, co2_rad_opt_ancil )
      co2_mix_ratio_now = 0.0_r_def
      co2_well_mixed = .false.
    end select


    ! Set MMR for Cs
    select case ( cs_rad_opt )
    case ( cs_rad_opt_off )
      cs_mix_ratio_now = 0.0_r_def
      cs_well_mixed = .true.
    case ( cs_rad_opt_constant )
      cs_mix_ratio_now = cs_mix_ratio
      cs_well_mixed = .true.
    case ( cs_rad_opt_time_varying )
      call gas_calc( cs_mix_ratio_now,    &
                     cs_clim_fcg_nyears,  &
                     cs_clim_fcg_years,   &
                     cs_clim_fcg_levls,   &
                     cs_clim_fcg_rates )
      cs_well_mixed = .true.
    case ( cs_rad_opt_prognostic, cs_rad_opt_ancil )
      cs_well_mixed = .false.
    end select

    ! Set MMR for H2
    select case ( h2_rad_opt )
    case ( h2_rad_opt_off )
      h2_mix_ratio_now = 0.0_r_def
      h2_well_mixed = .true.
    case ( h2_rad_opt_constant )
      h2_mix_ratio_now = h2_mix_ratio
      h2_well_mixed = .true.
    case ( h2_rad_opt_time_varying )
      call gas_calc( h2_mix_ratio_now,    &
                     h2_clim_fcg_nyears,  &
                     h2_clim_fcg_years,   &
                     h2_clim_fcg_levls,   &
                     h2_clim_fcg_rates )
      h2_well_mixed = .true.
    case ( h2_rad_opt_prognostic, h2_rad_opt_ancil )
      h2_well_mixed = .false.
    end select

    ! Set MMR for H2O
    select case ( h2o_rad_opt )
    case ( h2o_rad_opt_off )
      h2o_mix_ratio_now = 0.0_r_def
      h2o_well_mixed = .true.
    case ( h2o_rad_opt_constant )
      h2o_mix_ratio_now = h2o_mix_ratio
      h2o_well_mixed = .true.
    case ( h2o_rad_opt_time_varying )
      call gas_calc( h2o_mix_ratio_now,    &
                     h2o_clim_fcg_nyears,  &
                     h2o_clim_fcg_years,   &
                     h2o_clim_fcg_levls,   &
                     h2o_clim_fcg_rates )
      h2o_well_mixed = .true.
    case ( h2o_rad_opt_prognostic, h2o_rad_opt_ancil )
      h2o_well_mixed = .false.
    end select

    ! Set MMR for HCFC22
    select case ( hcfc22_rad_opt )
    case ( hcfc22_rad_opt_off )
      hcfc22_mix_ratio_now = 0.0_r_def
    case ( hcfc22_rad_opt_constant )
      hcfc22_mix_ratio_now = hcfc22_mix_ratio
    case ( hcfc22_rad_opt_time_varying )
      call gas_calc( hcfc22_mix_ratio_now,    &
                     hcfc22_clim_fcg_nyears,  &
                     hcfc22_clim_fcg_years,   &
                     hcfc22_clim_fcg_levls,   &
                     hcfc22_clim_fcg_rates )
    end select

    ! Set MMR for HCN
    select case ( hcn_rad_opt )
    case ( hcn_rad_opt_off )
      hcn_mix_ratio_now = 0.0_r_def
    case ( hcn_rad_opt_constant )
      hcn_mix_ratio_now = hcn_mix_ratio
    case ( hcn_rad_opt_time_varying )
      call gas_calc( hcn_mix_ratio_now,    &
                     hcn_clim_fcg_nyears,  &
                     hcn_clim_fcg_years,   &
                     hcn_clim_fcg_levls,   &
                     hcn_clim_fcg_rates )
      hcn_well_mixed = .true.
    case ( hcn_rad_opt_prognostic, hcn_rad_opt_ancil )
      hcn_well_mixed = .false.
    end select

    ! Set MMR for He
    select case ( he_rad_opt )
    case ( he_rad_opt_off )
      he_mix_ratio_now = 0.0_r_def
      he_well_mixed = .true.
    case ( he_rad_opt_constant )
      he_mix_ratio_now = he_mix_ratio
      he_well_mixed = .true.
    case ( he_rad_opt_time_varying )
      call gas_calc( he_mix_ratio_now,    &
                     he_clim_fcg_nyears,  &
                     he_clim_fcg_years,   &
                     he_clim_fcg_levls,   &
                     he_clim_fcg_rates )
      he_well_mixed = .true.
    case ( he_rad_opt_prognostic, he_rad_opt_ancil )
      he_well_mixed = .false.
    end select

    ! Set MMR for HFC134a
    select case ( hfc134a_rad_opt )
    case ( hfc134a_rad_opt_off )
      hfc134a_mix_ratio_now = 0.0_r_def
    case ( hfc134a_rad_opt_constant )
      hfc134a_mix_ratio_now = hfc134a_mix_ratio
    case ( hfc134a_rad_opt_time_varying )
      call gas_calc( hfc134a_mix_ratio_now,    &
                     hfc134a_clim_fcg_nyears,  &
                     hfc134a_clim_fcg_years,   &
                     hfc134a_clim_fcg_levls,   &
                     hfc134a_clim_fcg_rates )
    end select

    ! Set MMR for K
    select case ( k_rad_opt )
    case ( k_rad_opt_off )
      k_mix_ratio_now = 0.0_r_def
      k_well_mixed = .true.
    case ( k_rad_opt_constant )
      k_mix_ratio_now = k_mix_ratio
      k_well_mixed = .true.
    case ( k_rad_opt_time_varying )
      call gas_calc( k_mix_ratio_now,    &
                     k_clim_fcg_nyears,  &
                     k_clim_fcg_years,   &
                     k_clim_fcg_levls,   &
                     k_clim_fcg_rates )
      k_well_mixed = .true.
    case ( k_rad_opt_prognostic, k_rad_opt_ancil )
      k_well_mixed = .false.
    end select

    ! Set MMR for Li
    select case ( li_rad_opt )
    case ( li_rad_opt_off )
      li_mix_ratio_now = 0.0_r_def
      li_well_mixed = .true.
    case ( li_rad_opt_constant )
      li_mix_ratio_now = li_mix_ratio
      li_well_mixed = .true.
    case ( li_rad_opt_time_varying )
      call gas_calc( li_mix_ratio_now,    &
                     li_clim_fcg_nyears,  &
                     li_clim_fcg_years,   &
                     li_clim_fcg_levls,   &
                     li_clim_fcg_rates )
      li_well_mixed = .true.
    case ( li_rad_opt_prognostic, li_rad_opt_ancil )
      li_well_mixed = .false.
    end select

    ! Set MMR for N2
    select case ( n2_rad_opt )
    case ( n2_rad_opt_off )
      n2_mix_ratio_now = 0.0_r_def
      n2_well_mixed = .true.
    case ( n2_rad_opt_constant )
      n2_mix_ratio_now = n2_mix_ratio
      n2_well_mixed = .true.
    case ( n2_rad_opt_time_varying )
      call gas_calc( n2_mix_ratio_now,    &
                     n2_clim_fcg_nyears,  &
                     n2_clim_fcg_years,   &
                     n2_clim_fcg_levls,   &
                     n2_clim_fcg_rates )
      n2_well_mixed = .true.
    case ( n2_rad_opt_prognostic, n2_rad_opt_ancil )
      n2_well_mixed = .false.
    end select

    ! Set MMR for N2O
    select case ( n2o_rad_opt )
    case ( n2o_rad_opt_off )
      n2o_mix_ratio_now = 0.0_r_def
      n2o_well_mixed = .true.
    case ( n2o_rad_opt_constant )
      n2o_mix_ratio_now = n2o_mix_ratio
      n2o_well_mixed = .true.
    case ( n2o_rad_opt_time_varying )
      call gas_calc( n2o_mix_ratio_now,    &
                     n2o_clim_fcg_nyears,  &
                     n2o_clim_fcg_years,   &
                     n2o_clim_fcg_levls,   &
                     n2o_clim_fcg_rates )
      n2o_well_mixed = .true.
    case ( n2o_rad_opt_prognostic, n2o_rad_opt_ancil )
      n2o_well_mixed = .false.
    end select

    ! Set MMR for Na
    select case ( na_rad_opt )
    case ( na_rad_opt_off )
      na_mix_ratio_now = 0.0_r_def
      na_well_mixed = .true.
    case ( na_rad_opt_constant )
      na_mix_ratio_now = na_mix_ratio
      na_well_mixed = .true.
    case ( na_rad_opt_time_varying )
      call gas_calc( na_mix_ratio_now,    &
                     na_clim_fcg_nyears,  &
                     na_clim_fcg_years,   &
                     na_clim_fcg_levls,   &
                     na_clim_fcg_rates )
      na_well_mixed = .true.
    case ( na_rad_opt_prognostic, na_rad_opt_ancil )
      na_well_mixed = .false.
    end select

    ! Set MMR for NH3
    select case ( nh3_rad_opt )
    case ( nh3_rad_opt_off )
      nh3_mix_ratio_now = 0.0_r_def
      nh3_well_mixed = .true.
    case ( nh3_rad_opt_constant )
      nh3_mix_ratio_now = nh3_mix_ratio
      nh3_well_mixed = .true.
    case ( nh3_rad_opt_time_varying )
      call gas_calc( nh3_mix_ratio_now,    &
                     nh3_clim_fcg_nyears,  &
                     nh3_clim_fcg_years,   &
                     nh3_clim_fcg_levls,   &
                     nh3_clim_fcg_rates )
      nh3_well_mixed = .true.
    case ( nh3_rad_opt_prognostic, nh3_rad_opt_ancil )
      nh3_well_mixed = .false.
    end select

    ! Set MMR for O2
    select case ( o2_rad_opt )
    case ( o2_rad_opt_off )
      o2_mix_ratio_now = 0.0_r_def
      o2_well_mixed = .true.
    case ( o2_rad_opt_constant )
      o2_mix_ratio_now = o2_mix_ratio
      o2_well_mixed = .true.
    case ( o2_rad_opt_time_varying )
      call gas_calc( o2_mix_ratio_now,    &
                     o2_clim_fcg_nyears,  &
                     o2_clim_fcg_years,   &
                     o2_clim_fcg_levls,   &
                     o2_clim_fcg_rates )
      o2_well_mixed = .true.
    case ( o2_rad_opt_prognostic, o2_rad_opt_ancil )
      o2_well_mixed = .false.
    end select

    ! Set MMR for O3
    select case ( o3_rad_opt )
    case ( o3_rad_opt_off )
      o3_mix_ratio_now = 0.0_r_def
      o3_well_mixed = .true.
    case ( o3_rad_opt_constant )
      o3_mix_ratio_now = o3_mix_ratio
      o3_well_mixed = .true.
    case ( o3_rad_opt_time_varying )
      call gas_calc( o3_mix_ratio_now,    &
                     o3_clim_fcg_nyears,  &
                     o3_clim_fcg_years,   &
                     o3_clim_fcg_levls,   &
                     o3_clim_fcg_rates )
      o3_well_mixed = .true.
    case ( o3_rad_opt_prognostic, o3_rad_opt_profile, o3_rad_opt_ancil )
      o3_well_mixed = .false.
    end select

    ! Set MMR for Rb
    select case ( rb_rad_opt )
    case ( rb_rad_opt_off )
      rb_mix_ratio_now = 0.0_r_def
      rb_well_mixed = .true.
    case ( rb_rad_opt_constant )
      rb_mix_ratio_now = rb_mix_ratio
      rb_well_mixed = .true.
    case ( rb_rad_opt_time_varying )
      call gas_calc( rb_mix_ratio_now,    &
                     rb_clim_fcg_nyears,  &
                     rb_clim_fcg_years,   &
                     rb_clim_fcg_levls,   &
                     rb_clim_fcg_rates )
      rb_well_mixed = .true.
    case ( rb_rad_opt_prognostic, rb_rad_opt_ancil )
      rb_well_mixed = .false.
    end select

    ! Set MMR for SO2
    select case ( so2_rad_opt )
    case ( so2_rad_opt_off )
      so2_mix_ratio_now = 0.0_r_def
      so2_well_mixed = .true.
    case ( so2_rad_opt_constant )
      so2_mix_ratio_now = so2_mix_ratio
      so2_well_mixed = .true.
    case ( so2_rad_opt_time_varying )
      call gas_calc( so2_mix_ratio_now,    &
                     so2_clim_fcg_nyears,  &
                     so2_clim_fcg_years,   &
                     so2_clim_fcg_levls,   &
                     so2_clim_fcg_rates )
      so2_well_mixed = .true.
    case ( so2_rad_opt_prognostic, so2_rad_opt_ancil )
      so2_well_mixed = .false.
    end select

    ! Set MMR for TiO
    select case ( tio_rad_opt )
    case ( tio_rad_opt_off )
      tio_mix_ratio_now = 0.0_r_def
      tio_well_mixed = .true.
    case ( tio_rad_opt_constant )
      tio_mix_ratio_now = tio_mix_ratio
      tio_well_mixed = .true.
    case ( tio_rad_opt_time_varying )
      call gas_calc( tio_mix_ratio_now,    &
                     tio_clim_fcg_nyears,  &
                     tio_clim_fcg_years,   &
                     tio_clim_fcg_levls,   &
                     tio_clim_fcg_rates )
      tio_well_mixed = .true.
    case ( tio_rad_opt_prognostic, tio_rad_opt_ancil )
      tio_well_mixed = .false.
    end select

    ! Set MMR for VO
    select case ( vo_rad_opt )
    case ( vo_rad_opt_off )
      vo_mix_ratio_now = 0.0_r_def
      vo_well_mixed = .true.
    case ( vo_rad_opt_constant )
      vo_mix_ratio_now = vo_mix_ratio
      vo_well_mixed = .true.
    case ( vo_rad_opt_time_varying )
      call gas_calc( vo_mix_ratio_now,    &
                     vo_clim_fcg_nyears,  &
                     vo_clim_fcg_years,   &
                     vo_clim_fcg_levls,   &
                     vo_clim_fcg_rates )
      vo_well_mixed = .true.
    case ( vo_rad_opt_prognostic, vo_rad_opt_ancil )
      vo_well_mixed = .false.
    end select

  end subroutine gas_calc_all

end module gas_calc_all_mod
