!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of moisture fields
!!          at a given point based upon a specified analytic formula
module analytic_moisture_profiles_mod

use constants_mod,                only : r_def, i_def, pi
use log_mod,                      only : log_event,                &
                                         log_scratch_space,        &
                                         LOG_LEVEL_ERROR
use idealised_config_mod,         only : test_grabowski_clark,      &
                                         test_deep_baroclinic_wave, &
                                         test_isot_dry_atm
use physics_common_mod,           only : qsaturation
use planet_config_mod,            only : recip_epsilon, scaled_radius
use coord_transform_mod,          only : xyz2llr, central_angle
use base_mesh_config_mod,         only : geometry, &
                                         geometry_spherical
use initial_density_config_mod,   only : r1, x1, y1, z1, r2, x2, y2, z2
use deep_baroclinic_wave_mod,     only : deep_baroclinic_wave

implicit none

private

public :: analytic_moisture

contains

!> @brief Compute an analytic moisture field
!> @param[in] chi         Position in Cartesian coordinates
!> @param[in] temperature Air temperature in K
!> @param[in] pressure    Air pressure in Pa
!> @param[in] choice      Integer defining which specified formula to use
!> @result moisture The resulting moisture field
function analytic_moisture(chi, temperature, pressure, choice) result(moisture)

  implicit none
  real(kind=r_def),    intent(in) :: chi(3)
  real(kind=r_def),    intent(in) :: temperature
  real(kind=r_def),    intent(in) :: pressure
  integer(kind=i_def), intent(in) :: choice
  real(kind=r_def)                :: moisture
  real(kind=r_def)                :: r, r1, r2, xc, zc  ! Spatial distances
  real(kind=r_def)                :: h0, rel_hum        ! Relative humidities
  real(kind=r_def)                :: mr_sat             ! Saturation value
  real(kind=r_def)                :: long, lat, radius, l1, l2
  real(kind=r_def)                :: theta, rho, exner, u, v, w ! Dummy fields

  ! We'll always need mr_sat so find it right away. Pressure is needed in mbar
  mr_sat = qsaturation(temperature, 0.01_r_def*pressure)

  if ( geometry == geometry_spherical ) then
    call xyz2llr(chi(1), chi(2), chi(3), long, lat, radius)
    call central_angle(long, lat, x1, y1, l1)
    call central_angle(long, lat, x2, y2, l2)
  else
    long = chi(1)
    lat  = chi(2)
    l1 = sqrt((long-x1)**2 + (lat-y1)**2)
    l2 = sqrt((long-x2)**2 + (lat-y2)**2)
  end if

  select case( choice )

  case( test_deep_baroclinic_wave )
    call deep_baroclinic_wave(long, lat, radius-scaled_radius, &
                              exner, theta, rho,               &
                              u, v, w, moisture)
  ! Test from Grabowski and Clark (1991)
  ! Returns the relative humidity field
  case( test_grabowski_clark )
    ! Parameters
    xc = 0.0_r_def    ! central x position of bubble
    zc = 800.0_r_def  ! central z position of bubble
    r1 = 300.0_r_def  ! outer radius of bubble
    r2 = 200.0_r_def  ! inner radius of bubble
    h0 = 0.2_r_def    ! background relative humidity

    r = sqrt((chi(1)-xc)**2.0_r_def + (chi(3)-zc)**2.0_r_def)
    if ( r <= r2 ) then
      rel_hum = 1.0_r_def
    else if ( r <= r1 ) then
      rel_hum = h0 + (1.0_r_def - h0) &
                      *( cos( pi*(r-r2) / (2.0_r_def*(r1-r2)) ) )**2.0_r_def
    else
      rel_hum = h0
    end if

    ! Now invert relative humidity expression to get water vapour mixing ratio
    moisture = rel_hum * mr_sat / &
               ( 1.0_r_def + (1.0_r_def - rel_hum)*mr_sat*recip_epsilon)

  ! Isothermal **dry** atmosphere
  case( test_isot_dry_atm )
    moisture = 0.0_r_def

  case default
    ! In other cases, mixing ratio is set just under saturation value
    moisture = 0.99_r_def*mr_sat

  end select

end function analytic_moisture

end module analytic_moisture_profiles_mod
