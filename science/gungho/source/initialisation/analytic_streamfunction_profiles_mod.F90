!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_streamfunction_profiles_mod

use constants_mod,           only: r_def, i_def, pi
use extrusion_config_mod,    only: domain_top
use initial_wind_config_mod, only: profile_sbr_streamfunction,      &
                                   profile_dcmip301_streamfunction, &
                                   profile_div_free_reversible,     &
                                   profile_eternal_fountain,        &
                                   profile_rotational
use planet_config_mod,       only: scaled_radius
use log_mod,                 only: log_event,                &
                                   log_scratch_space,        &
                                   LOG_LEVEL_ERROR

implicit none

private
public :: analytic_streamfunction

contains

!> @brief Compute an analytic streamfunction field
!> @param[in] chi Position in physical coordinates
!> @param[in] choice Integer defining which specified formula to use
!> @param[in] num_options Number of scalar options to supply
!> @param[in] option Array of real values used to generate the initial profile
!> @param[in] time A real used for time-varying stream functions
!> @param[in] domain_max_x Domain maximum x-coordinte.
!> @return psi Result streamfunction field vector u = curl(psi)
function analytic_streamfunction( chi, choice, num_options,    &
                                  option, time, domain_max_x ) &
                          result( psi )

  implicit none

  real(kind=r_def),    intent(in) :: chi(3)
  integer(kind=i_def), intent(in) :: choice, num_options
  real(kind=r_def),    intent(in) :: option(num_options)
  real(kind=r_def),    intent(in) :: time
  real(kind=r_def),    intent(in) :: domain_max_x

  real(kind=r_def)                :: psi(3)
  real(kind=r_def)                :: s
  real(kind=r_def)                :: lat_pole, lon_pole
  real(kind=r_def)                :: u0, time_period, L
  real(kind=r_def)                :: vortex_zcentre, la, lb, ld, lr
  real(kind=r_def)                :: coeffs(4)


  select case ( choice )

    case ( profile_sbr_streamfunction, &
           profile_dcmip301_streamfunction)
      s = 0.5_r_def*(chi(3)/scaled_radius + 1.0_r_def)
      ! Turn off the height variation for the dcmip test
      if ( choice == profile_dcmip301_streamfunction ) s = 1.0_r_def

      lat_pole = pi/2.0_r_def - option(2)*pi
      lon_pole = pi/2.0_r_def + option(3)*pi

      psi(1) = 0.0_r_def
      psi(2) = 0.0_r_def
      psi(3) = s * option(1) * scaled_radius * &
               (sin(lat_pole)*sin(chi(2)) - cos(lat_pole)*sin(chi(1)-lon_pole)*cos(chi(2)) )

  case ( profile_eternal_fountain )
    ! The eternal fountain profile from Zerroukat & Allen 2020 SLIC paper
    u0 = 0.002_r_def * domain_max_x

    psi(1) = 0.0_r_def
    psi(3) = 0.0_r_def
    psi(2) = - u0 * domain_top / pi *                                             &
              sin(2.0_r_def * pi * (chi(1) / (2.0 * domain_max_x) + 0.5_r_def)) * &
              sin(pi * chi(3) / domain_top)
  case ( profile_div_free_reversible )
    ! A deformational divergence-free time-varying flow
    ! based upon a test in Cotter & Kuzmin 2016
    L = 2.0_r_def * domain_max_x
    time_period = L  ! Solution should be the same as initial condition at t = L

    psi(1) = 0.0_r_def
    psi(3) = 0.0_r_def
    psi(2) = L * domain_top / time_period *                                    &
              ( - chi(3) / domain_top                                          &
                + 1.0_r_def / pi * (-0.5_r_def + time / time_period)           &
                * sin(2.0_r_def * pi * (chi(1) / L - time / time_period))      &
                * sin(2.0_r_def * pi * chi(3) / domain_top))

  case ( profile_rotational )
    ! A solid body rotation in a vertical slice
    ! The stream function is smoothed towards the edge of the domain
    ! to prevent errors at the boundaries
    time_period = domain_top  ! One rotation when T = H
    vortex_zcentre = domain_top / 2.0_r_def
    lr = domain_top
    la = 10.0_r_def * lr / 25.0_r_def
    lb = 12.0_r_def * lr / 25.0_r_def
    ld = lr * sqrt((chi(1) / (domain_max_x * 2.0_r_def)) ** 2.0_r_def &
                   + ((chi(3) - vortex_zcentre) / domain_top) ** 2.0_r_def)

    coeffs(1) = pi / time_period
    coeffs(2) = - pi * la / (time_period * (lb - la))
    coeffs(3) = pi * la * lb / time_period
    coeffs(4) = pi * la * lb / time_period

    psi(1) = 0.0_r_def
    psi(3) = 0.0_r_def
    if (ld <= la) then
      psi(2) = coeffs(1) * ld ** 2.0_r_def
    else if (ld <= lb) then
      psi(2) = coeffs(2) * (ld - lb) ** 2.0_r_def + coeffs(3)
    else
      psi(2) = coeffs(4)
    end if

    case default
      write( log_scratch_space, '(A)' )  'Invalid streamfunction profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

end function analytic_streamfunction

end module analytic_streamfunction_profiles_mod
