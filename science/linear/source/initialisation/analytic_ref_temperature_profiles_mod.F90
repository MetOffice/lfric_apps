!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Functions to compute a variety of unperturbed analytic profiles.
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula.
module analytic_ref_temperature_profiles_mod

use constants_mod,         only: r_def, i_def
use log_mod,               only: log_event,                    &
                                 log_scratch_space,            &
                                 LOG_LEVEL_ERROR
use idealised_config_mod,  only: test_cold_bubble_x,           &
                                 test_cold_bubble_y,           &
                                 test_gaussian_hill,           &
                                 test_cosine_hill,             &
                                 test_slotted_cylinder,        &
                                 test_gravity_wave,            &
                                 test_warm_bubble,             &
                                 test_warm_bubble_3d,          &
                                 test_solid_body_rotation,     &
                                 test_solid_body_rotation_alt, &
                                 test_deep_baroclinic_wave,    &
                                 test_cos_phi,                 &
                                 test_cosine_bubble,           &
                                 test_div_free_reversible,     &
                                 test_eternal_fountain,        &
                                 test_curl_free_reversible,    &
                                 test_rotational,              &
                                 test_translational,           &
                                 test_vertical_cylinder

use reference_profile_mod, only: reference_profile

implicit none

private

public :: analytic_ref_temperature

contains

!> @brief Compute an analytic temperature field.
!> @param[in] chi    Position in physical coordinates
!> @param[in] choice Integer defining which specified formula to use
!> @return           The result temperature field
function analytic_ref_temperature(chi, choice) result(temperature)

  implicit none

  real(kind=r_def),    intent(in) :: chi(3)
  integer(kind=i_def), intent(in) :: choice
  real(kind=r_def)                :: temperature
  real(kind=r_def)                :: pressure, density

  temperature = 0.0_r_def
  call reference_profile(pressure, density, temperature, chi, choice)

end function analytic_ref_temperature

end module analytic_ref_temperature_profiles_mod
