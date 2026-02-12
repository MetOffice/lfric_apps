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

module kessler_kernel_mod

  use argument_mod,                  only: arg_type,                    &
                                           GH_FIELD, GH_WRITE, GH_READ, &
                                           CELL_COLUMN, GH_REAL, GH_SCALAR, &
                                            ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,                 only: r_def, i_def
  use driver_water_constants_mod,    only: Lv => latent_heat_h2o_condensation
  use fs_continuity_mod,             only: Wtheta
  use kernel_mod,                    only: kernel_type
  use physics_common_mod,            only: qsaturation
  use planet_config_mod,             only: recip_epsilon, kappa, cp, Rd, p_zero

  implicit none

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: kessler_kernel_type
      private
      type(arg_type) :: meta_args(9) = (/                                      &
          arg_type(GH_FIELD,   GH_REAL, GH_WRITE, WTHETA),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  WTHETA),                     &
          arg_type(GH_FIELD*6, GH_REAL, GH_WRITE, WTHETA),                     &
          arg_type(GH_FIELD*6, GH_REAL, GH_READ,  WTHETA),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  WTHETA),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  WTHETA),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  WTHETA),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
          arg_type(GH_SCALAR,  GH_REAL, GH_READ)                               &
          /)
      integer :: operates_on = CELL_COLUMN

  contains
      procedure, nopass :: kessler_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: kessler_code

contains

  !> @brief Performs a simple condensation/evaporation scheme with latent heating
  !! @param[in] nlayers Integer the number of layers
  !! @param[in,out] theta_inc     Potential temperature increment
  !! @param[in]     theta_n       Potential temperature in
  !! @param[in,out] mr_v_inc      Water vapour mixing ratio
  !! @param[in,out] mr_cl_inc     Liquid cloud mixing ratio
  !! @param[in,out] mr_r_inc      Rain mixing ratio
  !! @param[in,out] mr_s_inc      Snow mixing ratio
  !! @param[in,out] mr_g_inc      Graupel mixing ratio
  !! @param[in,out] mr_ci_inc     Ice cloud mixing ratio
  !! @param[in]     mr_v_n        Water vapour mixing ratio
  !! @param[in]     mr_cl_n       Liquid cloud mixing ratio
  !! @param[in]     mr_r_n        Rain mixing ratio
  !! @param[in]     mr_s_n        Snow mixing ratio
  !! @param[in]     mr_g_n        Graupel mixing ratio
  !! @param[in]     mr_ci_n       Ice cloud mixing ratio
  !! @param[in]     exner  Exner pressure at Wtheta points
  !! @param[in]     ndf_wtheta    Number of DoFs per cell for Wtheta
  !! @param[in]     undf_wtheta   Universal number of DoFs for wtheta
  !! @param[in]     map_wtheta    Integers mapping DoFs to columns for Wtheta
  subroutine kessler_code(nlayers, theta_inc, theta,             &
                          mr_v_inc, mr_cl_inc, mr_r_inc,         &
                          mr_s_inc, mr_g_inc, mr_ci_inc,         &
                          mr_v, mr_cl, mr_r,                     &
                          mr_s, mr_g, mr_ci, exner,       &
                          rho, height_at_wth, rain2d,    &
                          dt, ndf_wtheta, undf_wtheta, map_wtheta, &
                          ndf_w3_2d, undf_w3_2d, map_w3_2d)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta
    integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
    integer(kind=i_def), dimension(ndf_wtheta),  intent(in)    :: map_wtheta
    integer(kind=i_def), dimension(ndf_w3_2d),   intent(in)    :: map_w3_2d
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: theta_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: theta
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: exner
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: rho
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: height_at_wth
    real(kind=r_def),    dimension(undf_w3_2d),  intent(inout) :: rain2d
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: mr_v_inc, mr_cl_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: mr_r_inc, mr_ci_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: mr_s_inc, mr_g_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: mr_v, mr_cl
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: mr_r, mr_ci
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: mr_s, mr_g
    real(kind=r_def), intent(in) :: dt

    ! Internal variables
    integer(kind=i_def) :: df, wt_idx
    real(kind=r_def)    :: theta_np1(nlayers+1), mr_v_np1(nlayers+1)
    real(kind=r_def)    :: mr_cl_np1(nlayers+1), mr_r_np1(nlayers+1)
    real(kind=r_def)    :: temperature, pressure(nlayers+1)
    real(kind=r_def)    :: sed(nlayers+1), rain2d_acc
    real(kind=r_def)    :: mr_sat, dm_v, dm_r, Vr(nlayers+1), Rv

    real(kind=r_def)    :: Er, C, dt0
    real(kind=r_def)    :: time_counter
    real(kind=r_def), parameter :: k1 = 0.001_r_def
    real(kind=r_def), parameter :: k2 = 2.2_r_def
    real(kind=r_def), parameter :: a = 0.001_r_def
    real(kind=r_def), parameter :: rhoqr = 1000.0_r_def

    ! Set max number of iterations for converging scheme
    wt_idx = map_wtheta(1)

    Rv = Rd * recip_epsilon

    ! Initialise arrays
    pressure(:) = p_zero * exner(wt_idx : wt_idx+nlayers) ** (1.0_r_def/kappa)
    theta_np1(:) = theta(wt_idx : wt_idx+nlayers)
    mr_v_np1(:) = mr_v(wt_idx : wt_idx+nlayers)
    mr_cl_np1(:) = mr_cl(wt_idx : wt_idx+nlayers)
    mr_r_np1(:) = mr_r(wt_idx : wt_idx+nlayers)

    ! ======================================================================== !
    ! Calculate terminal velocity
    ! ======================================================================== !

    do df = 1, nlayers + 1
      Vr(df) = 36.34_r_def * SQRT(rho(wt_idx) / rho(wt_idx + df - 1))          &
        * (MAX(0.0_r_def, mr_r_np1(df)) * 0.001_r_def * rho(wt_idx + df - 1)) ** 0.1346_r_def
    end do

    ! ======================================================================== !
    ! Compute maximum time step
    ! ======================================================================== !

    dt0 = dt
    do df = 1, nlayers - 1
      ! NB: Original test for Vr /= 0 numerically unstable
      if (abs(Vr(df)) > 1.0E-12_r_def) then
         dt0 = min(dt0, 0.8_r_def*(height_at_wth(wt_idx+df) - height_at_wth(wt_idx+df-1)) / Vr(df))
      end if
    end do

    ! ======================================================================== !
    ! Loop through microphysics iteration
    ! ======================================================================== !

    ! time counter keeps track of the elapsed time during the subcycling process
    time_counter = 0.0_r_def

    ! initialize time-weighted accumulated precipitation
    rain2d_acc = 0.0_r_def
    ! Subcycle through the Kessler moisture processes,
    ! time loop ends when the physics time step is reached (within a margin of 1e-5 s)
    do while ( abs(dt - time_counter) > 1.0E-5_r_def )

      ! Precipitation rate (m_water/s) over the subcycled time step
      rain2d(map_w3_2d(1)) = rho(wt_idx) * mr_r_np1(wt_idx) * Vr(1) / rhoqr

      ! accumulate the preciptation rate over the subcycled time steps
      ! (weighted with the subcycled time step), unit is m_water
      rain2d_acc = rain2d_acc + rain2d(map_w3_2d(1)) * dt0

      ! Mass-weighted sedimentation term using upstream differencing
      do df = 1, nlayers
        sed(df) = dt0 *                                                        &
          ((rho(wt_idx+df) * mr_r_np1(df+1) * Vr(df+1))                        &
            - (rho(wt_idx+df-1) * mr_r_np1(df) * Vr(df)))                      &
          / (rho(wt_idx+df-1) * (height_at_wth(wt_idx+df) - height_at_wth(wt_idx+df-1)))
      end do
      sed(nlayers+1) = -dt0 * mr_r_np1(wt_idx+nlayers) * Vr(nlayers+1)         &
        / (0.5_r_def * (height_at_wth(wt_idx+nlayers)-height_at_wth(wt_idx+nlayers-1)))

      do df = 1, nlayers+1

        ! Autoconversion and collection rates following Klemp and Wilhelmson (1978), Eqs. (2.13a,b)
        ! the collection process is handled with a semi-implicit time stepping approach
        dm_r = mr_cl_np1(df)                                                   &
          - (mr_cl_np1(df) - dt0 * MAX(k1 * (mr_cl_np1(df) - a), 0.0_r_def))   &
          / (1.0_r_def + dt0 * k2 * mr_r_np1(df)**0.875_r_def)
        mr_cl_np1(df) = MAX(mr_cl_np1(df) - dm_r, 0.0_r_def)
        mr_r_np1(df) = MAX(mr_r_np1(df) + dm_r + sed(df), 0.0_r_def)

        ! Teten's formula
        temperature = theta_np1(df) * exner(wt_idx+df-1)
        ! This function takes pressure in mbar so divide by 100
        mr_sat = qsaturation(temperature, 0.01_r_def*pressure(df))

        ! Determine difference to saturation amount for vapour
        dm_v = (mr_v_np1(df) - mr_sat) /                                       &
                (1.0_r_def + (mr_sat * 4093.0_r_def * Lv / cp) /               &
                             (temperature - 36.0_r_def)**2)

        ! Evaporation of rain
        C = 1.6_r_Def + 124.9_r_def * (rho(wt_idx+df-1)*mr_r_np1(df))**0.2046_r_def
        Er = (1.0_r_def - mr_v_np1(df) / mr_sat)*C*(rho(wt_idx+df-1)*mr_r_np1(df))**0.525_r_def       &
            / (rho(wt_idx+df-1) * 5.4e5_r_def + 2.55e6_r_def / (pressure(df)*mr_sat))
        Er = MIN(dt0*Er, MAX(-dm_v - mr_cl_np1(df), mr_r_np1(df)))

        ! Clip to prevent negative cloud forming
        if (dm_v < 0.0_r_def) then
          dm_v = max(dm_v, -mr_cl_np1(df))
        end if

        ! Update fields
        mr_v_np1(df) = MAX(mr_v_np1(df) - dm_v + Er, 0.0_r_def)
        mr_cl_np1(df) = mr_cl_np1(df) + dm_v
        mr_r_np1(df) = MAX(mr_r_np1(df) - Er, 0.0_r_def)
        theta_np1(df) = theta_np1(df) + Lv / (cp * exner(wt_idx+df-1)) * (dm_v - Er)
      end do ! Loop over DoFs

      ! Compute the elapsed time
      time_counter = time_counter + dt0

      ! Recalculate liquid water terminal velocity (m/s)
      do df = 1, nlayers + 1
        Vr(df) = 36.34_r_def * SQRT(rho(wt_idx) / rho(wt_idx + df - 1))          &
          * (MAX(0.0_r_def, mr_r_np1(df)) * 0.001_r_def * rho(wt_idx + df - 1)) ** 0.1346_r_def
      end do

      ! recompute the time step
      dt0 = max(dt -  time_counter, 0.0_r_def)
      do df = 1, nlayers - 1
        ! NB: Original test for Vr /= 0 numerically unstable
        if (abs(Vr(df)) > 1.0E-12_r_def) then
          dt0 = min(dt0, 0.8_r_def*(height_at_wth(wt_idx+df) - height_at_wth(wt_idx+df-1)) / Vr(df))
        end if
      end do

    end do ! Microphysics time loop

  ! ========================================================================== !
  ! Calculate increments
  ! ========================================================================== !

  theta_inc(wt_idx : wt_idx + nlayers) = theta_np1(:) - theta(wt_idx : wt_idx + nlayers)
  mr_v_inc(wt_idx : wt_idx + nlayers) = mr_v_np1(:) - mr_v(wt_idx : wt_idx + nlayers)
  mr_cl_inc(wt_idx : wt_idx + nlayers) = mr_cl_np1(:) - mr_cl(wt_idx : wt_idx + nlayers)
  mr_r_inc(wt_idx : wt_idx + nlayers) = mr_r_np1(:) - mr_r(wt_idx : wt_idx + nlayers)

  end subroutine kessler_code

end module kessler_kernel_mod
