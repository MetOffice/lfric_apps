!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes the moisture source term for the potential temperature
!> @details In the presence of moisture, the potential temperature equation
!!          includes a source term proportional to the wind divergence:
!!          \f[
!!          \frac{\partial \theta}{\partial t} + \mathbf{u} \cdot \nabla \theta
!!          = \frac{1}{c_{vm}} \left(R_m - R_d \frac{c_{pm}}{c_{pd}}\right)
!!          \theta \left(\nabla \cdot \mathbf{u}\right)
!!          \f]
!!          which disappears in the absence of moisture, or if the specific heat
!!          capacities of water species are neglected.
!!          This kernel performs a pointwise computation of this source term at
!!          Wtheta points, using a Crank-Nicolson time discretisation to update
!!          the dry potential temperature field
module theta_moist_source_kernel_mod

  use argument_mod,               only : arg_type, GH_SCALAR, GH_FIELD,        &
                                         GH_REAL, GH_READWRITE, GH_READ, DOF
  use constants_mod,              only : r_def, EPS
  use driver_water_constants_mod, only : cpv => heat_capacity_h2o_vapour,      &
                                         Rv => gas_constant_h2o,               &
                                         cl => heat_capacity_h2o,              &
                                         ci => heat_capacity_ice
  use fs_continuity_mod,          only : Wtheta
  use kernel_mod,                 only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: theta_moist_source_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta),                    &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: theta_moist_source_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: theta_moist_source_code

contains

!> @brief Computes the moisture source term for the potential temperature
!> @param[in,out] theta        The dry potential temperature field
!> @param[in]     div_u        The wind divergence field
!> @param[in]     mr_v         Mixing ratio of water vapour
!> @param[in]     mr_cl        Mixing ratio of cloud liquid
!> @param[in]     mr_r         Mixing ratio of rain
!> @param[in]     mr_ci        Mixing ratio of cloud ice
!> @param[in]     mr_s         Mixing ratio of snow
!> @param[in]     mr_g         Mixing ratio of graupel
!> @param[in]     dt           Time step
!> @param[in]     cpd          Heat capacity of dry air at constant pressure
!> @param[in]     Rd           Gas constant for dry air
subroutine theta_moist_source_code( theta,       &
                                    div_u,       &
                                    mr_v,        &
                                    mr_cl,       &
                                    mr_r,        &
                                    mr_ci,       &
                                    mr_s,        &
                                    mr_g,        &
                                    dt,          &
                                    cpd,         &
                                    Rd )

  implicit none

  ! Arguments
  real(kind=r_def), intent(inout) :: theta
  real(kind=r_def), intent(in)    :: div_u
  real(kind=r_def), intent(in)    :: mr_v, mr_cl, mr_r
  real(kind=r_def), intent(in)    :: mr_s, mr_g, mr_ci
  real(kind=r_def), intent(in)    :: dt, cpd, Rd

  ! Internal variables
  real(kind=r_def) :: cvd, cvv, cpm, cvm, Rm, source_term

  ! Heat capacities and moist air gas constant
  Rm = Rd + mr_v * Rv
  cvd = cpd - Rd
  cvv = cpv - Rv
  cvm = cvd + mr_v * cvv + (mr_cl + mr_r) * cl + (mr_s + mr_g + mr_ci) * ci
  cpm = cvm + Rm

  ! Compute source term and update potential temperature
  source_term = div_u * (Rm - Rd * (cpm / cpd)) / cvm
  theta = theta * (1.0_r_def - 0.5_r_def * dt * source_term)                   &
                   / MAX(1.0_r_def + 0.5_r_def * dt * source_term, EPS)

end subroutine theta_moist_source_code

end module theta_moist_source_kernel_mod
