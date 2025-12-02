!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module for computing a linear hydrostatially balanced reference state
module reference_profile_mod

use base_mesh_config_mod,           only : geometry,           &
                                           geometry_spherical, &
                                           topology,           &
                                           topology_fully_periodic
use constants_mod,                  only : r_def, i_def
use coord_transform_mod,            only : xyz2llr
use generate_global_gw_fields_mod,  only : generate_global_gw_fields
use idealised_config_mod,           only : test_cold_bubble_x,    &
                                           test_cold_bubble_y,    &
                                           test_const_lapse_rate, &
                                           test_cosine_hill,      &
                                           test_gravity_wave,     &
                                           test_isentropic,       &
                                           test_isot_atm,         &
                                           test_isot_cold_atm,    &
                                           test_isot_dry_atm,     &
                                           test_warm_bubble,      &
                                           test_warm_bubble_3d,   &
                                           test_yz_cosine_hill
use initial_pressure_config_mod,    only : surface_pressure
use initial_temperature_config_mod, only : bvf_square, theta_surf
use planet_config_mod,              only : scaled_radius, gravity, Cp, Rd, &
                                           kappa, p_zero
use formulation_config_mod,         only : shallow
implicit none

private

public :: reference_profile

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the analytic reference profile at a single point
!! @param[out] exner_s   Pressure reference profile
!! @param[out] rho_s     Density reference profile
!! @param[out] theta_s   Potential temperature reference profile
!! @param[in] x          (x,y,z) coordinate field
!! @param[in] itest_option Choice of idealised profile
subroutine reference_profile(exner_s, rho_s, theta_s, x, itest_option)

implicit none

real(kind=r_def),    intent(in)           :: x(3)
integer(kind=i_def), intent(in)           :: itest_option
real(kind=r_def),    intent(out)          :: exner_s, rho_s, theta_s

real(kind=r_def), parameter :: lapse_rate     = 0.0065_r_def
real(kind=r_def)            :: nsq_over_g, z, u_s(3), lat, lon, r, t, p
real(kind=r_def)            :: exner_surf

exner_surf = (surface_pressure / p_zero) ** kappa

if ( geometry == geometry_spherical ) then  ! SPHERICAL DOMAIN
  call xyz2llr(x(1),x(2),x(3),lon,lat,r)
  z = r - scaled_radius

  select case( itest_option )

    case( test_const_lapse_rate )
      t = (theta_surf - lapse_rate * z)
      p = p_zero*(1.0_r_def - lapse_rate/theta_surf*z)**(gravity/(Rd*lapse_rate))

      ! Convert p (pressure) to exner_s (exner pressure)
      exner_s = (p/p_zero)**(kappa)

      ! Convert t (temperature) to theta (potential temperature)
      theta_s = t/exner_s

      rho_s = p/(Rd*t)

    case( test_isot_atm, test_isot_dry_atm )
      t = theta_surf
      p = p_zero * exp(-gravity * z / (Rd * t))

      ! Convert p (pressure) to exner_s (exner pressure)
      exner_s = (p / p_zero) ** (kappa)

      ! Convert t (temperature) to theta (potential temperature)
      theta_s = t / exner_s

      ! Calculate density using the ideal gas law
      rho_s = p / (Rd * t)

    case default
      ! Gravity wave test only for now
      call generate_global_gw_fields (lat, z, exner_s, u_s, theta_s, rho_s)
  end select

else                     ! PLANAR DOMAIN

  z = x(3)

  ! Calculate theta and exner for each biperiodic test
  select case( itest_option )
    case( test_gravity_wave, &
          test_isot_atm )
      nsq_over_g = bvf_square/gravity
      theta_s = theta_surf * exp ( nsq_over_g * z )
      exner_s = exner_surf - gravity**2/(Cp * theta_surf * bvf_square)   &
                   * (1.0_r_def - exp ( - nsq_over_g * z ))
    case( test_isot_cold_atm)
      theta_s = theta_surf * exp ( gravity / (theta_surf * cp) * z )
      exner_s = exner_surf * exp ( - gravity / (theta_surf * cp) * z )
    case( test_warm_bubble,    &
          test_cold_bubble_x,  &
          test_cold_bubble_y,  &   ! Density current test
          test_warm_bubble_3d, &
          test_isentropic )
      theta_s = theta_surf
      exner_s = exner_surf - gravity/(Cp*theta_surf)*z
    case( test_const_lapse_rate )
      theta_s = theta_surf * ((1.0_r_def - lapse_rate/theta_surf * z) &
                  **(1.0_r_def-gravity/(Cp*lapse_rate)))
      exner_s = exner_surf * ((1.0_r_def - lapse_rate/theta_surf * z) &
                  **(gravity/(Cp*lapse_rate)))
    !> @todo No values for the following idealised tests were provided and
    !>       this risked unexpected divide by zero errors. These errors are
    !>       avoided by setting to one. This keeps the trunk working but
    !>       these numbers have no scientific value.
    case (test_cosine_hill, &
          test_yz_cosine_hill)
      theta_s = 1.0_r_def
      exner_s = 1.0_r_def
    case default
      theta_s = theta_surf
      exner_s = exner_surf
  end select
  ! Calculate rho for all biperiodic tests
  rho_s   = p_zero/(Rd*theta_s) * exner_s ** ((1.0_r_def - kappa)/kappa)

end if

end subroutine reference_profile

end module reference_profile_mod
