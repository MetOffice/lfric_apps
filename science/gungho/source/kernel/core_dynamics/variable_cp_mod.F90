!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5.5 (Claude Code).
!-----------------------------------------------------------------------------
!> @brief Temperature-dependent specific heat capacity at constant pressure,
!!        and the potential temperature that is conserved under it.
!>
!> @details The heat capacity follows the power law
!!          \f[ c_p(T) = c_{p,ref} \left( \frac{T}{T_0} \right)^\nu \f]
!!          Integrating the dry adiabat under this law gives the potential
!!          temperature theta_cp_law, which is materially conserved in
!!          adiabatic motion:
!!          \f[ \theta_{cp}^\nu = T^\nu + \nu T_0^\nu
!!              \ln\left[ \left( \frac{p_0}{p} \right)^{R/c_{p,ref}} \right] \f]
!!          The exponent R/c_{p,ref} uses the law's own reference heat
!!          capacity. It is not the model's kappa, which is built on the
!!          namelist cp. theta_cp_law is also not the model's prognostic
!!          theta. With nu = 0 the law is a constant c_{p,ref}, and
!!          theta_cp_law reduces to the classical T (p_0/p)^{R/c_{p,ref}}.
!!          The law's parameters are passed in as arguments, so the functions
!!          hold no state and are safe to call from kernels.
module variable_cp_mod

  use constants_mod, only : r_def

  implicit none

  private

  public :: cp_law
  public :: theta_cp_law

contains

  !> @brief Heat capacity at constant pressure at a given temperature.
  !> @param[in] temperature      Temperature (K)
  !> @param[in] cp_law_ref       Heat capacity at the reference temperature
  !!                             (J/kg/K)
  !> @param[in] cp_law_t0        Reference temperature of the law (K)
  !> @param[in] cp_law_exponent  Exponent of the law
  !> @return    cp               Heat capacity at constant pressure (J/kg/K)
  elemental function cp_law( temperature, cp_law_ref, cp_law_t0, &
                             cp_law_exponent ) result( cp )

    implicit none

    real(kind=r_def), intent(in) :: temperature
    real(kind=r_def), intent(in) :: cp_law_ref, cp_law_t0, cp_law_exponent
    real(kind=r_def)             :: cp

    cp = cp_law_ref * ( temperature / cp_law_t0 )**cp_law_exponent

  end function cp_law

  !> @brief Potential temperature conserved along a dry adiabat under the
  !!        heat capacity law.
  !> @param[in] temperature      Temperature (K)
  !> @param[in] pressure         Pressure (Pa)
  !> @param[in] p_zero           Reference pressure (Pa)
  !> @param[in] rd               Specific gas constant (J/kg/K)
  !> @param[in] cp_law_ref       Heat capacity at the reference temperature
  !!                             (J/kg/K)
  !> @param[in] cp_law_t0        Reference temperature of the law (K)
  !> @param[in] cp_law_exponent  Exponent of the law
  !> @return    theta            Potential temperature (K)
  elemental function theta_cp_law( temperature, pressure, p_zero, rd,   &
                                   cp_law_ref, cp_law_t0, cp_law_exponent ) &
                                   result( theta )

    implicit none

    real(kind=r_def), intent(in) :: temperature, pressure, p_zero, rd
    real(kind=r_def), intent(in) :: cp_law_ref, cp_law_t0, cp_law_exponent
    real(kind=r_def)             :: theta

    real(kind=r_def) :: kappa_ref

    kappa_ref = rd / cp_law_ref

    if ( cp_law_exponent == 0.0_r_def ) then
      theta = temperature * ( p_zero / pressure )**kappa_ref
    else
      theta = ( temperature**cp_law_exponent                              &
                + cp_law_exponent * cp_law_t0**cp_law_exponent            &
                  * kappa_ref * log( p_zero / pressure ) )                &
              **( 1.0_r_def / cp_law_exponent )
    end if

  end function theta_cp_law

end module variable_cp_mod
