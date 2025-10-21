!-------------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief  Computes theta_m, the moist potential temperature
!> @details The moist potential temperature theta_m is materially conserved in
!!          the absence of phase changes. This kernel computes theta_m from the
!!          prognostic variables.
module dry_to_moist_theta_kernel_mod

  use argument_mod,               only : arg_type, GH_SCALAR, GH_FIELD,        &
                                         GH_REAL, GH_WRITE, GH_READ, DOF
  use constants_mod,              only : r_def
  use driver_water_constants_mod, only : cpv => heat_capacity_h2o_vapour,      &
                                         Rv => gas_constant_h2o,               &
                                         cl => heat_capacity_h2o,              &
                                         ci => heat_capacity_ice
  use fs_continuity_mod,          only : Wtheta
  use kernel_mod,                 only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: dry_to_moist_theta_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: dry_to_moist_theta_code
  end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: dry_to_moist_theta_code

contains

!> @brief Computes theta_m, the moist potential temperature, from theta_d
!> @param[in,out] theta_m    Moist potential temperature to compute
!> @param[in]     theta_d    Dry potential temperature
!> @param[in]     rho_at_wt  Density at Wtheta points
!> @param[in]     mr_v       Mixing ratio of water vapour
!> @param[in]     mr_cl      Mixing ratio of cloud liquid
!> @param[in]     mr_r       Mixing ratio of rain
!> @param[in]     mr_ci      Mixing ratio of cloud ice
!> @param[in]     mr_s       Mixing ratio of snow
!> @param[in]     mr_g       Mixing ratio of graupel
!> @param[in]     cpd        Heat capacity of dry air at constant pressure
!> @param[in]     Rd         Gas constant for dry air
!> @param[in]     p_zero     Reference pressure
subroutine dry_to_moist_theta_code( theta_m,     &
                                    theta_d,     &
                                    rho_at_wt,   &
                                    mr_v,        &
                                    mr_cl,       &
                                    mr_r,        &
                                    mr_ci,       &
                                    mr_s,        &
                                    mr_g,        &
                                    cpd,         &
                                    Rd,          &
                                    p_zero       )

  implicit none

  ! Arguments
  real(kind=r_def), intent(inout) :: theta_m
  real(kind=r_def), intent(in)    :: theta_d, rho_at_wt
  real(kind=r_def), intent(in)    :: mr_v, mr_cl, mr_r, mr_s, mr_g, mr_ci
  real(kind=r_def), intent(in)    :: cpd, Rd, p_zero

  ! Internal variables
  real(kind=r_def)    :: cvd, cvv, cpm, cvm, Rm, kappa_d, kappa_m, p, T

  ! Heat capacities and moist air gas constant
  Rm = Rd + mr_v * Rv
  cvd = cpd - Rd
  cvv = cpv - Rv
  cvm = cvd + mr_v * cvv + (mr_cl + mr_r) * cl + (mr_s + mr_g + mr_ci) * ci
  cpm = cvm + Rm

  ! Exponents
  kappa_d = Rd / cpd
  kappa_m = Rm / cpm

  ! Compute moist potential temperature
  T = (theta_d * (Rm * rho_at_wt / p_zero)**kappa_d)**(1.0_r_def / (1.0_r_def - kappa_d))
  p = rho_at_wt * Rm * T
  theta_m = T * (p / p_zero)**(-kappa_m)

end subroutine dry_to_moist_theta_code

end module dry_to_moist_theta_kernel_mod
