!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic geopotential profiles.
!> @details Collection of functions to return the value of the geopotential or
!!          the surface geopotential at a given point based upon a specified
!!          analytic formula for the shallow water miniapp.
module analytic_geopot_profiles_mod

  ! Constants
  use constants_mod,          only: i_def, r_def, pi
  use log_mod,                only: log_event,                &
                                    log_scratch_space,        &
                                    LOG_LEVEL_ERROR

  ! Functions
  use coord_transform_mod,    only: xyz2ll
  use galewsky_test_case_mod, only: galewsky_profile

  ! Configurations
  use base_mesh_config_mod,   only: geometry, &
                                    geometry_spherical
  use planet_config_mod,      only: scaled_radius, gravity, &
                                    scaled_omega, scaling_factor
  use shallow_water_settings_config_mod,                           &
                              only: ref_gp,                        &
                                    swe_test_swe_geostr_balance,   &
                                    swe_test_swe_geostr_imbalance, &
                                    swe_test_swe_gaussian_hill,    &
                                    swe_test_swe_mountain,         &
                                    swe_test_swe_galewsky,         &
                                    swe_test_swe_vortex_field,     &
                                    swe_test_swe_thermal_dbl_vortex
  use shallow_water_test_coeff_config_mod, &
                              only: x1, y1, xr, yr, mag1

  implicit none

contains

  !============================================================================
  !> @brief Compute an analytic geopotential field for the shallow water miniapp.
  !> @param[in] chi      Position in physical coordinates
  !> @param[in] choice   Which test to set up
  !> @param[in] domain_x Domain size in x-direction.
  !> @result    geopot The result geopotential field
  function analytic_geopot(chi, choice, domain_x) result(geopot)

    implicit none

    real(kind=r_def), intent(in) :: chi(3)
    integer,          intent(in) :: choice
    real(kind=r_def), intent(in) :: domain_x

    real(kind=r_def)             :: geopot
    real(kind=r_def)             :: long, lat
    real(kind=r_def)             :: fc, ubya, href, psi_dummy
    real(kind=r_def)             :: x1_vortex, x2_vortex, y1_vortex, y2_vortex, &
                                    sx, sy, ox, oy,                             &
                                    xc1, xc2, yc1, yc2, lx
    real(kind=r_def)             :: u0, gh0
    real(kind=r_def), parameter  :: day = 86400.0_r_def

    select case( choice )

    case ( swe_test_swe_geostr_balance )
      if ( geometry == geometry_spherical ) then
        call xyz2ll(chi(1),chi(2),chi(3),long,lat)
        u0 = 2.0_r_def * pi *scaled_radius / ( 12.0_r_def * day )
        gh0 = ref_gp
        geopot = scaling_factor*gh0 &
                - (scaled_radius*scaled_omega*u0 + u0**2_i_def / 2.0_r_def )*(sin(lat))**2_i_def
      else
        fc = 2.0_r_def * scaled_omega
        geopot = ref_gp + fc * cos ( 4.0_r_def * pi * ( chi(2) + 0.5_r_def ) ) / ( 4.0_r_def * pi )
      end if
    case (swe_test_swe_vortex_field)
        geopot = ref_gp
    case( swe_test_swe_gaussian_hill )
      if ( geometry == geometry_spherical ) then
        call xyz2ll(chi(1),chi(2),chi(3),long,lat)
        geopot = ref_gp + exp( - ( ( long - x1 ) / xr ) ** 2 &
                               - ( ( lat  - y1 ) / yr ) ** 2 ) * mag1
      else
        geopot = ref_gp + exp( - ( ( chi(1) - x1 ) / xr ) ** 2 &
                               - ( ( chi(2) - y1 ) / yr ) ** 2 ) * mag1
      end if
    case( swe_test_swe_geostr_imbalance )
      fc = 2.0_r_def * scaled_omega
      geopot = ref_gp + fc* sin ( 4.0_r_def * pi * ( chi(2) + 0.5_r_def ) ) / ( 4.0_r_def * pi )

    case( swe_test_swe_mountain )
      ubya = 20.0_r_def / scaled_radius
      href = 5960.0_r_def
      geopot = scaling_factor * href * gravity &
               - chi(3) ** 2 * ( scaled_omega * ubya + ubya ** 2 / 2.0_r_def )

    case ( swe_test_swe_galewsky )
      call xyz2ll(chi(1),chi(2),chi(3),long,lat)
      call galewsky_profile(lat, long, psi_dummy, geopot)

    case ( swe_test_swe_thermal_dbl_vortex )
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
      x1_vortex = lx * sin ( pi * ( chi(1) - xc1 ) / lx ) / ( pi * sx )
      x2_vortex = lx * sin ( pi * ( chi(1) - xc2 ) / lx ) / ( pi * sx )
      y1_vortex = lx * sin ( pi * ( chi(2) - yc1 ) / lx ) / ( pi * sy )
      y2_vortex = lx * sin ( pi * ( chi(2) - yc2 ) / lx ) / ( pi * sy )
      geopot = ref_gp - 0.1_r_def * ref_gp *                                   &
               (  exp( - 0.5_r_def * ( x1_vortex ** 2 + y1_vortex ** 2 ) )     &
                + exp( - 0.5_r_def * ( x2_vortex ** 2 + y2_vortex ** 2 ) )     &
                - 4.0_r_def * pi * sx * sy / ( lx ** 2 ) )

    case default
      write( log_scratch_space, '(A)' )  'Invalid geopotential profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

  end function analytic_geopot

  !============================================================================
  !> @brief Compute an analytic surface geopotential field for the shallow water miniapp.
  !> @param[in] chi      Position in physical coordinates
  !> @result    s_geopot The result surcace geopotential field
  function analytic_surface_geopot(chi, choice) result(s_geopot)

    implicit none

    real(kind=r_def), intent(in) :: chi(3)
    integer,          intent(in) :: choice
    real(kind=r_def)             :: s_geopot
    real(kind=r_def)             :: long, lat
    real(kind=r_def)             :: longc, latc
    real(kind=r_def)             :: rad, rad_ext, mheight

    select case( choice )

    case ( swe_test_swe_geostr_balance, swe_test_swe_gaussian_hill, swe_test_swe_galewsky, &
           swe_test_swe_geostr_imbalance, swe_test_swe_thermal_dbl_vortex, swe_test_swe_vortex_field )
      s_geopot = 0.0_r_def

    case( swe_test_swe_mountain )
      call xyz2ll(chi(1),chi(2),chi(3),long,lat)
      longc   = 1.0_r_def * pi / 2.0_r_def
      latc    = pi / 6.0_r_def
      rad_ext = pi / 9.0_r_def
      rad     = min( rad_ext, sqrt( ( long - longc ) ** 2 + ( lat - latc ) ** 2 ) )
      mheight = 2000.0_r_def
      s_geopot = mheight * gravity * (1.0_r_def - rad/rad_ext)

    case default
      write( log_scratch_space, '(A)' )  'Invalid surface geopotential profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

  end function analytic_surface_geopot

  !============================================================================
  !> @brief Compute an analytic tracer field for the shallow water miniapp.
  !> @param[in] chi    Position in physical coordinates
  !> @param[in] domain_x Domain size in x-direction.
  !> @result    tracer The result tracer field
  function analytic_tracer(chi, domain_x) result(tracer)

    implicit none

    real(kind=r_def), intent(in) :: chi(3)
    real(kind=r_def), intent(in) :: domain_x

    real(kind=r_def)             :: tracer
    real(kind=r_def)             :: long, lat

    if ( geometry == geometry_spherical ) then
      call xyz2ll(chi(1),chi(2),chi(3),long,lat)
      if ( abs(lat-pi/4.0_r_def) < pi/16.0_r_def ) then
        tracer = 1.0_r_def
      else
        tracer = 0.0_r_def
      end if
    else
      if ( abs(chi(2)) < domain_x/16.0_r_def ) then
        tracer = 1.0_r_def
      else
        tracer = 0.0_r_def
      end if
    end if

  end function analytic_tracer

end module analytic_geopot_profiles_mod
