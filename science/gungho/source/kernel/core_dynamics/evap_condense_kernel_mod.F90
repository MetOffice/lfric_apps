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

  use argument_mod,                  only: arg_type,                    &
                                           GH_FIELD, GH_WRITE, GH_READ, &
                                           CELL_COLUMN, GH_REAL
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
  type, public, extends(kernel_type) :: evap_condense_kernel_type
      private
      type(arg_type) :: meta_args(5) = (/                  &
          arg_type(GH_FIELD,   GH_REAL, GH_WRITE, WTHETA), &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  WTHETA), &
          arg_type(GH_FIELD*6, GH_REAL, GH_WRITE, WTHETA), &
          arg_type(GH_FIELD*6, GH_REAL, GH_READ,  WTHETA), &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  WTHETA)  &
          /)
      integer :: operates_on = CELL_COLUMN

  contains
      procedure, nopass :: evap_condense_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: evap_condense_code

contains

  !> @brief Performs a simple condensation/evaporation scheme with latent heating
  !! @param[in] nlayers Integer the number of layers
  !! @param[in,out] theta_inc     Potential temperature increment
  !! @param[in]     theta_n       Potential temperature in
  !! @param[in,out] mr_v_inc      Water vapour mixing ratio
  !! @param[in,out] mr_cl_inc     Liquid cloud mixing ratio
  !! @param[in,out] mr_r_inc      Rain mixing ratio
  !! @param[in,out] mr_ci_inc     Ice cloud mixing ratio
  !! @param[in,out] mr_s_inc      Snow mixing ratio
  !! @param[in,out] mr_g_inc      Graupel mixing ratio
  !! @param[in]     mr_v_n        Water vapour mixing ratio
  !! @param[in]     mr_cl_n       Liquid cloud mixing ratio
  !! @param[in]     mr_r_n        Rain mixing ratio
  !! @param[in]     mr_ci_n       Ice cloud mixing ratio
  !! @param[in]     mr_s_n        Snow mixing ratio
  !! @param[in]     mr_g_n        Graupel mixing ratio
  !! @param[in]     exner_at_wth  Exner pressure at Wtheta points
  !! @param[in]     ndf_wtheta    Number of DoFs per cell for Wtheta
  !! @param[in]     undf_wtheta   Universal number of DoFs for wtheta
  !! @param[in]     map_wtheta    Integers mapping DoFs to columns for Wtheta
  subroutine evap_condense_code(nlayers, theta_inc, theta_n,           &
                                mr_v_inc, mr_cl_inc, mr_r_inc,         &
                                mr_ci_inc, mr_s_inc, mr_g_inc,         &
                                mr_v_n, mr_cl_n, mr_r_n,               &
                                mr_ci_n, mr_s_n, mr_g_n, exner_at_wth, &
                                ndf_wtheta, undf_wtheta, map_wtheta    )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, undf_wtheta
    integer(kind=i_def), dimension(ndf_wtheta),  intent(in)    :: map_wtheta
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: theta_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: theta_n
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: exner_at_wth
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: mr_v_inc, mr_cl_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: mr_r_inc, mr_ci_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(inout) :: mr_s_inc, mr_g_inc
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: mr_v_n, mr_cl_n
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: mr_r_n, mr_ci_n
    real(kind=r_def),    dimension(undf_wtheta), intent(in)    :: mr_s_n, mr_g_n

    ! Internal variables
    integer(kind=i_def) :: k, df, lower_df, i, max_iter
    real(kind=r_def)    :: theta, mr_v, mr_cl
    real(kind=r_def)    :: mr_sat, dm_v, Rv
    real(kind=r_def)    :: temperature, pressure

    ! Set max number of iterations for converging scheme
    max_iter = 1

    Rv = Rd * recip_epsilon

    do k = 0, nlayers-1
      ! If we are at bottom do lower df, otherwise only top df
      if ( k == 0 ) then
        lower_df = 1
      else
        lower_df = 2
      end if

      do df = lower_df, 2

        ! First guesses
        theta = theta_n(map_wtheta(df)+k)
        mr_v = mr_v_n(map_wtheta(df)+k)
        mr_cl = mr_cl_n(map_wtheta(df)+k)
        pressure = p_zero * exner_at_wth(map_wtheta(df)+k) ** (1.0_r_def/kappa)

        ! Iterate towards a solution
        do i = 1, max_iter

          temperature = theta * exner_at_wth(map_wtheta(df)+k)
          ! This function takes pressure in mbar so divide by 100
          mr_sat = qsaturation(temperature, 0.01_r_def*pressure)

          ! Determine difference to saturation amount for vapour
          dm_v = (mr_v - mr_sat) /                                &
                 (1.0_r_def + (mr_sat * Lv ** 2.0_r_def) /        &
                              (cp * Rv * temperature ** 2.0_r_def))

          ! Clip to prevent negative cloud forming
          if (dm_v < 0.0_r_def) then
            dm_v = max(dm_v, -mr_cl_n(map_wtheta(df)+k))
          end if

          ! Update fields
          mr_v = mr_v - dm_v
          mr_cl = mr_cl + dm_v
          theta = theta * (1.0_r_def + dm_v * Lv / (cp * temperature))

        end do ! Thermodynamics iteration

        theta_inc(map_wtheta(df)+k) = theta - theta_n(map_wtheta(df)+k)
        mr_v_inc(map_wtheta(df)+k) = mr_v - mr_v_n(map_wtheta(df)+k)
        mr_cl_inc(map_wtheta(df)+k) = mr_cl - mr_cl_n(map_wtheta(df)+k)

      end do ! Loop over DoFs

    end do ! Loop over layers

  end subroutine evap_condense_code

end module evap_condense_kernel_mod
