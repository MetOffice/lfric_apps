module gas_calc_all_mod

  use constants_mod,       only: i_def, r_def, r_um
  use gas_calc_mod,        only: gas_calc
  use rad_input_mod,       only: co2_mmr
  use well_mixed_gases_config_mod, only: clim_fcg_years_ch4,      &
                                         clim_fcg_nyears_ch4,     &
                                         clim_fcg_levls_ch4,      &
                                         clim_fcg_rates_ch4,      &
                                         clim_fcg_years_co2,      &
                                         clim_fcg_nyears_co2,     &
                                         clim_fcg_levls_co2,      &
                                         clim_fcg_rates_co2,      &
                                         clim_fcg_years_n2o,      &
                                         clim_fcg_nyears_n2o,     &
                                         clim_fcg_levls_n2o,      &
                                         clim_fcg_rates_n2o,      &
                                         clim_fcg_years_cfc11,    &
                                         clim_fcg_nyears_cfc11,   &
                                         clim_fcg_levls_cfc11,    &
                                         clim_fcg_rates_cfc11,    &
                                         clim_fcg_years_cfc12,    &
                                         clim_fcg_nyears_cfc12,   &
                                         clim_fcg_levls_cfc12,    &
                                         clim_fcg_rates_cfc12,    &
                                         clim_fcg_years_cfc113,   &
                                         clim_fcg_nyears_cfc113,  &
                                         clim_fcg_levls_cfc113,   &
                                         clim_fcg_rates_cfc113,   &
                                         clim_fcg_years_hcfc22,   &
                                         clim_fcg_nyears_hcfc22,  &
                                         clim_fcg_levls_hcfc22,   &
                                         clim_fcg_rates_hcfc22,   &
                                         clim_fcg_years_hfc134a,  &
                                         clim_fcg_nyears_hfc134a, &
                                         clim_fcg_levls_hfc134a,  &
                                         clim_fcg_rates_hfc134a,  &
                                         cfc113_mix_ratio,        &
                                         cfc11_mix_ratio,         &
                                         cfc12_mix_ratio,         &
                                         ch4_mix_ratio,           &
                                         co2_mix_ratio,           &
                                         hcfc22_mix_ratio,        &
                                         hfc134a_mix_ratio,       &
                                         n2_mix_ratio,            &
                                         n2o_mix_ratio,           &
                                         so2_mix_ratio,           &
                                         l_time_varying_co2,      &
                                         l_time_varying_ch4,      &
                                         l_time_varying_n2o,      &
                                         l_time_varying_cfc11,    &
                                         l_time_varying_cfc12,    &
                                         l_time_varying_hcfc22,   &
                                         l_time_varying_hfc134a,  &
                                         l_time_varying_cfc113

  implicit none

  ! latest updated gas mixing ratios
  real(r_def), public :: co2_mix_ratio_now
  real(r_def), public :: n2o_mix_ratio_now
  real(r_def), public :: ch4_mix_ratio_now
  real(r_def), public :: cfc11_mix_ratio_now
  real(r_def), public :: cfc12_mix_ratio_now
  real(r_def), public :: cfc113_mix_ratio_now
  real(r_def), public :: hcfc22_mix_ratio_now
  real(r_def), public :: hfc134a_mix_ratio_now

contains

  subroutine gas_calc_all()

    implicit none

    ! Update C02
    if ( l_time_varying_co2 ) then
      call gas_calc( co2_mix_ratio_now,    &
                     clim_fcg_nyears_co2,  &
                     clim_fcg_years_co2,   &
                     clim_fcg_levls_co2,   &
                     clim_fcg_rates_co2 )
    else
      co2_mix_ratio_now = co2_mix_ratio
    end if
    ! CO2 value needed by JULES - contained stored in rad_input
    co2_mmr = real(co2_mix_ratio_now, r_um)

    ! Update N2O
    if ( l_time_varying_n2o ) then
      call gas_calc( n2o_mix_ratio_now,    &
                     clim_fcg_nyears_n2o,  &
                     clim_fcg_years_n2o,   &
                     clim_fcg_levls_n2o,   &
                     clim_fcg_rates_n2o )
    else
      n2o_mix_ratio_now = n2o_mix_ratio
    end if

    ! Update CH4
    if ( l_time_varying_ch4 ) then
      call gas_calc( ch4_mix_ratio_now,    &
                     clim_fcg_nyears_ch4,  &
                     clim_fcg_years_ch4,   &
                     clim_fcg_levls_ch4,   &
                     clim_fcg_rates_ch4 )
    else
      ch4_mix_ratio_now = ch4_mix_ratio
    end if

    ! Update CFC11
    if ( l_time_varying_cfc11 ) then
      call gas_calc( cfc11_mix_ratio_now,    &
                     clim_fcg_nyears_cfc11,  &
                     clim_fcg_years_cfc11,   &
                     clim_fcg_levls_cfc11,   &
                     clim_fcg_rates_cfc11 )
    else
      cfc11_mix_ratio_now = cfc11_mix_ratio
    end if

    ! Update CFC12
    if ( l_time_varying_cfc12 ) then
      call gas_calc( cfc12_mix_ratio_now,    &
                     clim_fcg_nyears_cfc12,  &
                     clim_fcg_years_cfc12,   &
                     clim_fcg_levls_cfc12,   &
                     clim_fcg_rates_cfc12 )
    else
      cfc12_mix_ratio_now = cfc12_mix_ratio
    end if

    ! Update CFC113
    if ( l_time_varying_cfc113 ) then
      call gas_calc( cfc113_mix_ratio_now,    &
                     clim_fcg_nyears_cfc113,  &
                     clim_fcg_years_cfc113,   &
                     clim_fcg_levls_cfc113,   &
                     clim_fcg_rates_cfc113 )
    else
      cfc113_mix_ratio_now = cfc113_mix_ratio
    end if

    ! Update HCFC22
    if ( l_time_varying_hcfc22 ) then
      call gas_calc( hcfc22_mix_ratio_now,    &
                     clim_fcg_nyears_hcfc22,  &
                     clim_fcg_years_hcfc22,   &
                     clim_fcg_levls_hcfc22,   &
                     clim_fcg_rates_hcfc22 )
    else
      hcfc22_mix_ratio_now = hcfc22_mix_ratio
    end if

    ! Update HFC134A
    if ( l_time_varying_hfc134a ) then
      call gas_calc( hfc134a_mix_ratio_now,    &
                     clim_fcg_nyears_hfc134a,  &
                     clim_fcg_years_hfc134a,   &
                     clim_fcg_levls_hfc134a,   &
                     clim_fcg_rates_hfc134a )
    else
      hfc134a_mix_ratio_now = hfc134a_mix_ratio
    end if

  end subroutine gas_calc_all

end module gas_calc_all_mod
