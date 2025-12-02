!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Computes theta_e, the wet equivalent potential temperature.
!>
!> @details Computes theta_e, the wet equivalent potential temperature, at
!>          Wtheta points from the potential temperature, the Exner pressure
!>          and the mixing ratio of water vapour.
!>
module theta_e_kernel_mod

  use argument_mod,               only : arg_type, GH_SCALAR, &
                                         GH_FIELD, GH_REAL,   &
                                         GH_WRITE, GH_READ,   &
                                         DOF
  use constants_mod,              only : r_def, i_def
  use driver_water_constants_mod, only : latent_heat_h2o_condensation
  use planet_config_mod,          only : Rd
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
  type, public, extends(kernel_type) :: theta_e_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: theta_e_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: theta_e_code

contains

!> @brief Computes theta_e, the wet equivalent potential temperature.
!> @param[in,out] theta_e        Wet equivalent potential temperature field
!> @param[in]     theta          Input potential temperature field
!> @param[in]     exner_at_wt    Exner pressure at Wtheta points
!> @param[in]     mr_v           Mixing ratio of water vapour
!> @param[in]     mr_cl          Mixing ratio of cloud liquid
!> @param[in]     cpd            Heat capacity of dry air at constant pressure
!> @param[in]     cpv            Heat capacity of water vap at constant pressure
!> @param[in]     cl             Heat capacity of liquid water
!> @param[in]     Rd             Gas constant for dry air
!> @param[in]     p_zero         Reference pressure
!> @param[in]     recip_epsilon  Reciprocal of the ratio of the gas constants
subroutine theta_e_code( theta_e,     &
                         theta,       &
                         exner_at_wt, &
                         mr_v,        &
                         mr_cl,       &
                         cpd,         &
                         cpv,         &
                         cl,          &
                         Rd,          &
                         p_zero,      &
                         recip_epsilon )

  implicit none

  ! Arguments
  real(kind=r_def), intent(inout) :: theta_e
  real(kind=r_def), intent(in)    :: theta
  real(kind=r_def), intent(in)    :: exner_at_wt
  real(kind=r_def), intent(in)    :: mr_v, mr_cl
  real(kind=r_def), intent(in)    :: cpd, cpv, cl, Rd, p_zero, recip_epsilon

  ! Internal variables
  real(kind=r_def) :: temperature, pressure, dry_pressure, kappa, Lv
  real(kind=r_def), parameter :: ref_temperature = 273.15_r_def

  temperature = theta * exner_at_wt
  pressure = p_zero * exner_at_wt**(cpd / Rd)
  dry_pressure = pressure / (1.0_r_def + recip_epsilon * mr_v)
  Lv = latent_heat_h2o_condensation - (cl - cpv)*(temperature - ref_temperature)
  kappa = Rd / (cpd + cl * (mr_v + mr_cl))

  theta_e = (                                                                  &
    temperature * (dry_pressure / p_zero) ** (-kappa)                          &
    * EXP( (Lv * mr_v) / (temperature * (cpd + cl * (mr_v + mr_cl))) )         &
  )

end subroutine theta_e_code

end module theta_e_kernel_mod
