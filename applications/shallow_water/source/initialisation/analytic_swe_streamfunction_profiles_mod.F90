!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic streamfunction profiles.
!> @details Collection of functions to return the value of a streamfunction at
!!          a given point based upon a specified analytic formula for the
!!          shallow water miniapp.
module analytic_swe_streamfunction_profiles_mod

  ! Constants
  use constants_mod,          only: r_def, pi
  use log_mod,                only: log_event,         &
                                    log_scratch_space, &
                                    LOG_LEVEL_ERROR

  ! Configuration
  use base_mesh_config_mod,   only: geometry, &
                                    geometry_spherical
  use planet_config_mod,      only: scaled_radius, scaled_omega
  use shallow_water_settings_config_mod,                   &
                              only: swe_test_swe_galewsky, &
                                    swe_test_swe_geostr_balance

  ! Functions
  use galewsky_test_case_mod, only: galewsky_profile

  implicit none

contains

  !> @brief Compute an analytic streamfunction field for the shallow water miniapp.
  !> @param[in] chi Position in physical coordinates
  !> @result    psi The result streamfunction field
  function analytic_swe_streamfunction(chi, choice) result(psi)

    implicit none
    real(kind=r_def), intent(in) :: chi(3)
    integer,          intent(in) :: choice
    real(kind=r_def)             :: psi(3)
    real(kind=r_def)             :: geopot_dummy
    real(kind=r_def)             :: u0
    real(kind=r_def), parameter  :: day = 86400.0_r_def

    psi(1:2) = 0.0_r_def

    select case( choice )

    case ( swe_test_swe_galewsky )
      ! Psi kernel passes cartesian or spherical coordinates depending on domain
      ! Arange input and psi such that the wind initialisation works correctly
      call galewsky_profile(chi(2), chi(1), psi(3), geopot_dummy)

    case ( swe_test_swe_geostr_balance )

      if ( geometry == geometry_spherical ) then

        u0 = 2.0_r_def * pi *scaled_radius / ( 12.0_r_def * day )
        psi(3) =  u0*scaled_radius*sin(chi(2))

      else

        psi(3) = -( cos(4.0_r_def*pi*(chi(2)+0.5_r_def)) )/(4.0_r_def*pi)

      end if

    case default
      write( log_scratch_space, '(A)' )  'Invalid streamfunction profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

  end function analytic_swe_streamfunction

end module analytic_swe_streamfunction_profiles_mod
