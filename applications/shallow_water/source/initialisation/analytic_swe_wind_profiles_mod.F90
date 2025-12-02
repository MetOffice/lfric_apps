!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic wind profiles.
!> @details Collection of functions to return the value of the wind field at a given
!!          point based upon a specified analytic formula for the shallow
!!          water miniapp.
module analytic_swe_wind_profiles_mod

  ! Constants
  use constants_mod,          only: r_def, pi
  use log_mod,                only: log_event,         &
                                    log_scratch_space, &
                                    LOG_LEVEL_ERROR

  ! Configuration
  use base_mesh_config_mod,   only: geometry, &
                                    geometry_spherical
  use planet_config_mod,      only: scaled_radius, scaled_omega
  use shallow_water_settings_config_mod,                           &
                              only: ref_gp,                        &
                                    swe_test_swe_geostr_balance,   &
                                    swe_test_swe_geostr_imbalance, &
                                    swe_test_swe_gaussian_hill,    &
                                    swe_test_swe_mountain,         &
                                    swe_test_swe_vortex_field,     &
                                    swe_test_swe_thermal_dbl_vortex
  use shallow_water_test_coeff_config_mod, &
                              only: u1, u2

  implicit none

contains

  !> @brief Compute an analytic velocity field for the shallow water miniapp.
  !> @param[in] chi  Position in physical coordinates
  !> @param[in] choice
  !> @param[in] domain_x Domain size in x-direction.
  !> @result    wind The resulting velocity field
  function analytic_swe_wind( chi, choice, domain_x ) result(wind)

    implicit none

    real(kind=r_def), intent(in)   :: chi(3)
    integer,          intent(in)   :: choice
    real(kind=r_def), intent(in)   :: domain_x
    real(kind=r_def), dimension(3) :: wind

    real(kind=r_def), parameter    :: day = 86400.0_r_def
    real(kind=r_def)               :: fc, lx, x1, x2, y1, y2, sx, sy, ox, oy, &
                                      x11, x22, y11, y22, xc1, xc2, yc1, yc2
    real(kind=r_def)               :: u0

    ! For spherical geometry chi from kernel is (lon, lat, rad)
    ! For planar geometry chi from kernel is (x, y, z)

    select case( choice )

    case ( swe_test_swe_geostr_balance )
      if ( geometry == geometry_spherical ) then
        ! Test 2 from Williamson et al. (JCP 1992)
        u0 = 2.0_r_def * pi *scaled_radius / ( 12.0_r_def * day )
        wind(1) = u0 * cos(chi(2))
        wind(2) = 0.0_r_def
        wind(3) = 0.0_r_def
      else
        wind(1) = sin ( 4.0_r_def * pi * ( chi(2) + 0.5_r_def ) )
        wind(2) = 0.0_r_def
        wind(3) = 0.0_r_def
      end if

    case ( swe_test_swe_vortex_field)
      if ( geometry == geometry_spherical ) then
        ! Test from Kent, Jablonowski, Thuburn and Wood (QJRMS 2016)
        wind(1) =  3.0_r_def/4.0_r_def*exp(-(chi(2)**8.0_r_def)/2.0_r_def)*( &
                   32.0_r_def*sin(4.0_r_def*chi(1))*cos(8.0_r_def*chi(2)) &
                   - 19.2_r_def*cos(3.0_r_def*chi(1))*sin(6.0_r_def*chi(2)) &
                   - 6.0_r_def*3.0_r_def*cos(5.0_r_def*chi(1))*sin(10.0_r_def*chi(2)) &
                   + cos(2.0_r_def*chi(2))  )
        wind(2) = -3.0_r_def/4.0_r_def*exp(-(chi(2)**8.0_r_def)/2.0_r_def)*( &
                  16.0_r_def*sin(8.0_r_def*chi(2))*cos(4.0_r_def*chi(1)) &
                  -9.6_r_def*cos(6.0_r_def*chi(2))*sin(3.0_r_def*chi(1)) &
                  -3.0_r_def*sin(5.0_r_def*chi(1))*cos(10.0_r_def*chi(2)) &
                  +0.5_r_def*cos(2.0_r_def*chi(1)) )
        wind(3) = 0.0_r_def
      else
        ! Test from Kent, Thuburn and Wood (QJRMS 2012)
        u0 = 10.0_r_def
        wind(1) = u0*sin(8.0_r_def*pi*chi(1))*cos(8.0_r_def*pi*chi(2))/(16.0_r_def*pi) &
                -u0*cos(6.0_r_def*pi*chi(1))*sin(6.0_r_def*pi*chi(2))/(30.0_r_def*pi) &
                -u0*3.0_r_def*cos(10.0_r_def*pi*chi(1))*sin(4.0_r_def*pi*chi(2))/(290.0_r_def*pi) &
                + u0*cos(2.0_r_def*pi*chi(2))/(100.0_r_def*pi)
        wind(2) = -u0*sin(8.0_r_def*pi*chi(2))*cos(8.0_r_def*pi*chi(1))/(16.0_r_def*pi) &
                +u0*cos(6.0_r_def*pi*chi(2))*sin(6.0_r_def*pi*chi(1))/(30.0_r_def*pi) &
                +u0*3.0_r_def*sin(10.0_r_def*pi*chi(1))*cos(4.0_r_def*pi*chi(2))/(116.0_r_def*pi) &
                - u0*cos(2.0_r_def*pi*chi(1))/(100.0_r_def*pi)
        wind(3) = 0.0_r_def
      end if

    case( swe_test_swe_gaussian_hill )
      if ( geometry == geometry_spherical ) then
        wind(1) = 0.0_r_def
        wind(2) = 0.0_r_def
        wind(3) = 0.0_r_def
      else
        wind(1) = u1
        wind(2) = u2
        wind(3) = 0.0_r_def
      end if

    case( swe_test_swe_geostr_imbalance )
      wind(1) = 0.0_r_def
      wind(2) = sin ( 2.0_r_def * pi * ( chi(1) + 0.5_r_def ) )
      wind(3) = 0.0_r_def

    case( swe_test_swe_mountain )
      ! Test 5 from Williamson et al. (JCP 1992) on spherical domain
      u0 = 20.0_r_def
      wind(1) = u0 * cos(chi(2))
      wind(2) = 0.0_r_def
      wind(3) = 0.0_r_def

    case ( swe_test_swe_thermal_dbl_vortex )
      fc  = 2.0_r_def * scaled_omega
      ! Assuming square planar domain
      lx = domain_x
      sx = lx * 3.0_r_def / 40.0_r_def
      sy = lx * 3.0_r_def / 40.0_r_def
      ox = 0.1_r_def
      oy = 0.1_r_def
      xc1 = - lx * ox
      xc2 =   lx * ox
      yc1 = - lx * oy
      yc2 =   lx * oy
      x1 = lx * sin ( pi * ( chi(1) - xc1 ) / lx ) / ( pi * sx )
      x2 = lx * sin ( pi * ( chi(1) - xc2 ) / lx ) / ( pi * sx )
      y1 = lx * sin ( pi * ( chi(2) - yc1 ) / lx ) / ( pi * sy )
      y2 = lx * sin ( pi * ( chi(2) - yc2 ) / lx ) / ( pi * sy )
      x11 = lx * sin ( 2.0_r_def * pi * ( chi(1) - xc1 ) / lx ) / ( pi * sx * 2.0_r_def )
      x22 = lx * sin ( 2.0_r_def * pi * ( chi(1) - xc2 ) / lx ) / ( pi * sx * 2.0_r_def  )
      y11 = lx * sin ( 2.0_r_def * pi * ( chi(2) - yc1 ) / lx ) / ( pi * sy * 2.0_r_def  )
      y22 = lx * sin ( 2.0_r_def * pi * ( chi(2) - yc2 ) / lx ) / ( pi * sy * 2.0_r_def  )

      wind(1) = - 0.1_r_def * ref_gp * (   y11 * exp( - 0.5_r_def * ( x1 ** 2 + y1 ** 2 ) ) &
                                         + y22 * exp( - 0.5_r_def * ( x2 ** 2 + y2 ** 2 ) ) ) / ( fc * sy )
      wind(2) =   0.1_r_def * ref_gp * (   x11 * exp( - 0.5_r_def * ( x1 ** 2 + y1 ** 2 ) ) &
                                         + x22 * exp( - 0.5_r_def * ( x2 ** 2 + y2 ** 2 ) ) ) / ( fc * sx )
      wind(3) =   0.0_r_def

    case default
      write( log_scratch_space, '(A)' )  'Invalid wind profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

end function analytic_swe_wind

end module analytic_swe_wind_profiles_mod
