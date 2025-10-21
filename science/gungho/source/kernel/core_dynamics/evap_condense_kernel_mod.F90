!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs a simple condensation/evaporation scheme with latent heating

!> @details Given the atmospheric temperature and pressure, this kernel computes
!!          the saturation mixing ratio of water vapour. Any excess vapour is
!!          condensed to cloud liquid, while any cloud liquid in an unsaturated
!!          environment is evaporated to water vapour. The potential temperature
!!          is adjusted to capture the effects of the latent heat release or
!!          absorption associated with this phase change.
!!          Note: this only works with the lowest order spaces

module evap_condense_kernel_mod

  use argument_mod,                only: arg_type, GH_SCALAR,         &
                                         GH_FIELD, GH_WRITE, GH_READ, &
                                         DOF, GH_REAL
  use constants_mod,               only: r_def, i_def
  use driver_water_constants_mod,  only: Lv0 => latent_heat_h2o_condensation
  use fs_continuity_mod,           only: Wtheta
  use kernel_mod,                  only: kernel_type
  use physics_common_mod,          only: qsaturation
  use planet_config_mod,           only: recip_epsilon, kappa, cpd => cp, Rd, p_zero

  implicit none

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: evap_condense_kernel_type
    private
    type(arg_type) :: meta_args(13) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, WTHETA),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                        &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = DOF

  contains
      procedure, nopass :: evap_condense_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: evap_condense_code

contains

  !> @brief Performs a simple condensation/evaporation scheme with latent heating
  !> @param[in,out] theta_inc    Potential temperature increment
  !> @param[in]     theta_n      Potential temperature input
  !> @param[in,out] mr_v_inc     Water vapour mixing ratio increment
  !> @param[in,out] mr_cl_inc    Liquid cloud mixing ratio increment
  !> @param[in]     mr_v_n       Water vapour mixing ratio input
  !> @param[in]     mr_cl_n      Liquid cloud mixing ratio input
  !> @param[in]     exner_at_wt  Exner pressure at Wtheta points
  !> @param[in]     Rd           Gas constant for dry air
  !> @param[in]     Rv           Gas constant for water vapour
  !> @param[in]     cpd          Heat capacity of dry air at constant pressure
  !> @param[in]     cpv          Heat capacity of water vap at constant pressure
  !> @param[in]     cl           Heat capacity of liquid water
  !> @param[in]     p_zero       Reference pressure
  subroutine evap_condense_code( theta_inc, theta_n,                           &
                                 mr_v_inc, mr_cl_inc,                          &
                                 mr_v_n, mr_cl_n,                              &
                                 exner_at_wt,                                  &
                                 Rd, Rv, cpd, cpv, cl, p_zero )

    implicit none

    ! Arguments
    real(kind=r_def), intent(inout) :: theta_inc
    real(kind=r_def), intent(in)    :: theta_n
    real(kind=r_def), intent(in)    :: exner_at_wt
    real(kind=r_def), intent(inout) :: mr_v_inc, mr_cl_inc
    real(kind=r_def), intent(in)    :: mr_v_n, mr_cl_n
    real(kind=r_def), intent(in)    :: Rd, Rv, cpd, cpv, cl, p_zero

    ! Internal variables
    real(kind=r_def) :: theta_np1, mr_v_np1, mr_cl_np1
    real(kind=r_def) :: cvd, cvv, kappa
    real(kind=r_def) :: mr_sat, dm_v, Lv, Rm, cpm, cvm
    real(kind=r_def) :: temperature, pressure
    real(kind=r_def), parameter :: ref_temperature = 273.15_r_def

    ! Convert to temperature and pressure
    kappa = cpd / Rd
    pressure = p_zero * exner_at_wt ** (1.0_r_def/kappa)
    temperature = theta_n * exner_at_wt

    ! Thermodynamic quantities
    kappa = Rd / cpd
    cvd = cpd - Rd
    cvv = cpv - Rv
    cpm = cpd + mr_v_n * cpv + mr_cl_n * cl
    cvm = cvd + mr_v_n * cvv + mr_cl_n * cl
    Rm = Rd + mr_v_n * Rv
    Lv = Lv0 - (cl - cpv)*(temperature - ref_temperature)

    ! This function takes pressure in mbar so divide by 100
    mr_sat = qsaturation(temperature, 0.01_r_def*pressure)

    ! Calculate vapour increment, which will be used for other increments
    dm_v = (mr_v_n - mr_sat) /                                                 &
            (1.0_r_def + (mr_sat * Lv ** 2.0_r_def) /                          &
                          (cpd * Rv * temperature ** 2.0_r_def))

    ! Clip to prevent negative cloud forming
    if (dm_v < 0.0_r_def) then
      dm_v = MAX(dm_v, -mr_cl_n)
    end if

    ! Update fields
    mr_v_np1 = mr_v_n - dm_v
    mr_cl_np1 = mr_cl_n + dm_v
    theta_np1 = theta_n * (                                                    &
        1.0_r_def + dm_v * ((cvd * Lv / (cvm * cpd * temperature))             &
        - (Rv / cvm) * (1 - ((Rd * cpm) / (cpd * Rm))))                        &
    )

    ! Compute final increments
    theta_inc = theta_np1 - theta_n
    mr_v_inc = - dm_v
    mr_cl_inc = mr_cl_np1 - mr_cl_n

  end subroutine evap_condense_code

end module evap_condense_kernel_mod
