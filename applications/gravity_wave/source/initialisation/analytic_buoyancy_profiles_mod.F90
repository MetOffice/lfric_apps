!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_buoyancy_profiles_mod

use constants_mod,                only : r_def, pi
use log_mod,                      only : log_event,                &
                                         log_scratch_space,        &
                                         LOG_LEVEL_ERROR
use coord_transform_mod,           only : xyz2llr
use base_mesh_config_mod,          only : geometry,           &
                                          geometry_spherical
use planet_config_mod,             only : scaled_radius, gravity
use generate_global_gw_fields_mod, only : generate_global_gw_pert
use reference_profile_mod,         only : reference_profile
use idealised_config_mod,          only : test_gravity_wave

implicit none

contains

!> @brief Compute an analytic buoyancy field
!> @param[in] chi Position in physical coordinates
!> @result buoyancy The result buoyancy field
function analytic_buoyancy(chi) result(buoyancy)

  implicit none
  real(kind=r_def), intent(in) :: chi(3)
  real(kind=r_def)             :: buoyancy

  real(kind=r_def), parameter  :: b0 = 0.01_r_def
  real(kind=r_def), parameter  :: XC     = -15000.0_r_def
  real(kind=r_def), parameter  :: A      = 5000.0_r_def
  real(kind=r_def), parameter  :: H      = 10000.0_r_def
  real(kind=r_def)             :: long, lat, radius
  real(kind=r_def)             :: pressure, density
  real(kind=r_def)             :: theta_0, theta_p

  if ( geometry == geometry_spherical ) then
    call xyz2llr(chi(1),chi(2),chi(3),long,lat,radius)
    call reference_profile(pressure, density, theta_0, chi, test_gravity_wave)
    theta_p = generate_global_gw_pert(long,lat,radius-scaled_radius)

    buoyancy = theta_p/theta_0 * gravity
  else
    buoyancy = b0 * sin ( pi * chi(3) / H ) &
              / ( 1.0_r_def + ( chi(1) - XC )**2/A**2 )
  end if

end function analytic_buoyancy

end module analytic_buoyancy_profiles_mod
