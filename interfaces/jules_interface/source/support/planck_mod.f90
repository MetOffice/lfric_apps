!----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Module with function to calculate the planck function.

module planck_mod
  implicit none

  contains

  function planck(t, lambda)

  use constants_mod,                   only : r_def

  use science_chemistry_constants_mod, only : boltzmann,       &
                                              planck_constant, &
                                              speed_of_light

  implicit none

  ! Name of function: planck (unit: W m^-3 steradian^-1)
  real(r_def) :: planck
  ! Temperature (unit: Kelvin)
  real(r_def), intent(in) :: t
  ! Wavelength (unit: Metre)
  real(r_def), intent(in) :: lambda

  ! Local variables.
  ! Inverse exponential factor
  real(r_def) :: exponential

  ! Evaluate a negative exponential to ensure conditioning.
  exponential = exp( - (planck_constant * speed_of_light / boltzmann) / (lambda * t) )

  ! planck function
  planck = 2.0_r_def * planck_constant * speed_of_light**2 * exponential &
  / ( (lambda**5) * (1.0_r_def - exponential) )

  end

end module planck_mod
