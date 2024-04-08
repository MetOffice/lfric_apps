!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic buoyancy profiles.
!> @details Collection of functions to return the value of the buoyancy field at
!!          a given point based upon a specified analytic formula for the
!!          shallow water miniapp.
module analytic_swe_buoyancy_profiles_mod

  ! Constants
  use constants_mod,          only: r_def, pi
  use log_mod,                only: log_event,         &
                                    log_scratch_space, &
                                    LOG_LEVEL_ERROR

  ! Configuration
  use shallow_water_settings_config_mod,                           &
                              only: ref_gp, thermal_swe,           &
                                    swe_test_swe_geostr_balance,   &
                                    swe_test_swe_geostr_imbalance, &
                                    swe_test_swe_gaussian_hill,    &
                                    swe_test_swe_mountain,         &
                                    swe_test_swe_galewsky,         &
                                    swe_test_swe_vortex_field,     &
                                    swe_test_swe_thermal_dbl_vortex

  implicit none

  contains

  !> @brief Compute an analytic buoyancy field for the shallow water miniapp.
  !> @param[in] chi      Position in physical coordinates
  !> @param[in] choice
  !> @param[in] domain_x Domain size in x-direction.
  !> @result    buoyancy The result buoyancy field
  function analytic_swe_buoyancy(chi, choice, domain_x) result(buoyancy)

    use analytic_geopot_profiles_mod, only: analytic_geopot

    implicit none

    real(kind=r_def), intent(in) :: chi(3)
    integer,          intent(in) :: choice
    real(kind=r_def), intent(in) :: domain_x
    real(kind=r_def)             :: buoyancy
    real(kind=r_def)             :: lx, gp

    select case( choice )

    case ( swe_test_swe_geostr_imbalance, swe_test_swe_mountain, swe_test_swe_galewsky )
      ! Constant 1 buoyancy recovers the shallow water equations
      buoyancy = 1.0_r_def

    case ( swe_test_swe_vortex_field )
      if ( thermal_swe ) then
        gp = analytic_geopot(chi, swe_test_swe_gaussian_hill, domain_x)
        buoyancy = ( gp / ref_gp ) ** 2
      else
        buoyancy = 1.0_r_def
      end if

    case ( swe_test_swe_thermal_dbl_vortex )
      ! Assuming square planar domain
      lx = domain_x
      buoyancy = 1.0_r_def + 0.05_r_def * sin( 2.0_r_def * pi * chi(1) / lx )

    case ( swe_test_swe_geostr_balance )
      if ( thermal_swe ) then
        gp = analytic_geopot(chi, swe_test_swe_geostr_balance, domain_x)
        buoyancy = 1.0_r_def + 0.05_r_def * ( ref_gp / gp ) ** 2
      else
        buoyancy = 1.0_r_def
      end if

    case ( swe_test_swe_gaussian_hill )
      if ( thermal_swe ) then
        gp = analytic_geopot(chi, swe_test_swe_gaussian_hill, domain_x)
        buoyancy = ( gp / ref_gp ) ** 2
      else
        buoyancy = 1.0_r_def
      end if

    case default
      write( log_scratch_space, '(A)' )  'Invalid buoyancy profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

  end function analytic_swe_buoyancy

end module analytic_swe_buoyancy_profiles_mod
