!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5.5 (Claude Code).
!-----------------------------------------------------------------------------
!> @brief Computes the source term for the potential temperature that arises
!!        from a temperature-dependent heat capacity
!> @details Potential temperature and the Exner pressure are defined with the
!!          constant heat capacity c_{p0}. When the true heat capacity c_p(T)
!!          differs from it, the dry potential temperature equation gains a
!!          source term proportional to the wind divergence:
!!          \f[
!!          \frac{\partial \theta}{\partial t} + \mathbf{u} \cdot \nabla \theta
!!          = \frac{R \left(c_p(T) - c_{p0}\right)}{c_{p0} \left(c_p(T) - R\right)}
!!          \theta \left(\nabla \cdot \mathbf{u}\right)
!!          \f]
!!          which vanishes where c_p(T) = c_{p0}. Where c_p(T) > c_{p0},
!!          expansion raises theta; where c_p(T) < c_{p0}, it lowers theta.
!!          c_p(T) follows the power law in variable_cp_mod. The temperature
!!          comes from theta and the density through the equation of state.
!!          This kernel performs a pointwise computation of this source term at
!!          Wtheta points, using a Crank-Nicolson time discretisation to update
!!          the dry potential temperature field.
module theta_variable_cp_source_kernel_mod

  use argument_mod,      only : arg_type, GH_SCALAR, GH_FIELD,                 &
                                GH_REAL, GH_READWRITE, GH_READ, DOF
  use constants_mod,     only : r_def, EPS
  use fs_continuity_mod, only : Wtheta
  use kernel_mod,        only : kernel_type
  use variable_cp_mod,   only : cp_law

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: theta_variable_cp_source_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: theta_variable_cp_source_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: theta_variable_cp_source_code

contains

!> @brief Computes the variable heat capacity source term for the potential
!!        temperature
!> @param[in,out] theta      The dry potential temperature field
!> @param[in]     div_u      The wind divergence field
!> @param[in]     rho_at_wt  The dry density at Wtheta points
!> @param[in]     dt         The time step
!> @param[in]     cp         The constant heat capacity on which theta and
!!                           exner are defined
!> @param[in]     Rd         The gas constant
!> @param[in]     p_zero     The reference pressure
!> @param[in]     cp_law_ref The heat capacity of the c_p(T) law at its
!!                           reference temperature
!> @param[in]     cp_law_t0  The reference temperature of the c_p(T) law
!> @param[in]     cp_law_exponent The exponent of the c_p(T) law
subroutine theta_variable_cp_source_code( theta,       &
                                          div_u,       &
                                          rho_at_wt,   &
                                          dt,          &
                                          cp,          &
                                          Rd,          &
                                          p_zero,      &
                                          cp_law_ref,  &
                                          cp_law_t0,   &
                                          cp_law_exponent )

  implicit none

  ! Arguments
  real(kind=r_def), intent(inout) :: theta
  real(kind=r_def), intent(in)    :: div_u
  real(kind=r_def), intent(in)    :: rho_at_wt
  real(kind=r_def), intent(in)    :: dt, cp, Rd, p_zero
  real(kind=r_def), intent(in)    :: cp_law_ref, cp_law_t0, cp_law_exponent

  ! Internal variables
  real(kind=r_def) :: kappa, temperature, cp_of_t, source_term

  ! Temperature from the equation of state
  kappa = Rd / cp
  temperature = (theta * (Rd * rho_at_wt / p_zero)**kappa)                     &
                **(1.0_r_def / (1.0_r_def - kappa))

  cp_of_t = cp_law( temperature, cp_law_ref, cp_law_t0, cp_law_exponent )

  ! Compute source term and update potential temperature.
  ! source_term is the coefficient S in D(theta)/Dt = -S theta
  source_term = - div_u * Rd * (cp_of_t - cp) / (cp * (cp_of_t - Rd))
  theta = theta * (1.0_r_def - 0.5_r_def * dt * source_term)                   &
                   / MAX(1.0_r_def + 0.5_r_def * dt * source_term, EPS)

end subroutine theta_variable_cp_source_code

end module theta_variable_cp_source_kernel_mod
