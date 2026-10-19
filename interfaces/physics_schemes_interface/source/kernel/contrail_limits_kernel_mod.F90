!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code (Claude Opus 5.5), 2026-10-15.

!> @brief Calculate the pressures of the lower and upper contrail limits.

module contrail_limits_kernel_mod

  use argument_mod,      only: arg_type,                  &
                               GH_FIELD, GH_SCALAR,       &
                               GH_READ, GH_WRITE,         &
                               GH_REAL, CELL_COLUMN,      &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def, rmdi
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: contrail_limits_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                     &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: contrail_limits_code
  end type contrail_limits_kernel_type

  public :: contrail_limits_code

contains

  !> @brief   Find the pressures of the lower and upper contrail limits.
  !> @details Contrails are expected where the air is colder than the
  !>          Mintra -14 degC environment curve, T = a*p**b, an empirical
  !>          temperature-pressure threshold for contrail formation behind
  !>          jet aircraft. The lower limit is the surface if contrails
  !>          form at the lowest level, otherwise the first crossing of the
  !>          curve; the upper limit is the last crossing above it. Only
  !>          crossings above 50 hPa pressure and below 18 km height are
  !>          searched. Each crossing pressure is found assuming a constant
  !>          lapse rate across the layer. Limits not found are set to rmdi.
  !> @param[in]     nlayers   The number of layers
  !> @param[in]     theta     Potential temperature (K)
  !> @param[in]     exner     Exner pressure on Wtheta
  !> @param[in]     height    Height above planet_radius on Wtheta (m)
  !> @param[in,out] p_lower   Pressure of the lower contrail limit (Pa)
  !> @param[in,out] p_upper   Pressure of the upper contrail limit (Pa)
  !> @param[in]     p_zero    Reference surface pressure (Pa)
  !> @param[in]     kappa     Ratio of Rd to cp
  !> @param[in]     g_over_r  Gravity divided by the dry-air gas constant
  !> @param[in]     ndf_wth   Number of degrees of freedom per cell (Wtheta)
  !> @param[in]     undf_wth  Number of unique degrees of freedom (Wtheta)
  !> @param[in]     map_wth   Dofmap for the cell at the base of the column
  !> @param[in]     ndf_2d    Number of degrees of freedom per cell (2D)
  !> @param[in]     undf_2d   Number of unique degrees of freedom (2D)
  !> @param[in]     map_2d    Dofmap for the 2D cell
  subroutine contrail_limits_code(nlayers,                    &
                                  theta, exner, height,       &
                                  p_lower, p_upper,           &
                                  p_zero, kappa, g_over_r,    &
                                  ndf_wth, undf_wth, map_wth, &
                                  ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d

    real(kind=r_def), intent(in),    dimension(undf_wth) :: theta
    real(kind=r_def), intent(in),    dimension(undf_wth) :: exner
    real(kind=r_def), intent(in),    dimension(undf_wth) :: height
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: p_lower
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: p_upper
    real(kind=r_def), intent(in) :: p_zero
    real(kind=r_def), intent(in) :: kappa
    real(kind=r_def), intent(in) :: g_over_r

    ! Mintra -14 degC curve, T = mintra_coeff * p**mintra_power (T in K,
    ! p in Pa)
    real(kind=r_def), parameter :: mintra_coeff = 137.14816_r_def
    real(kind=r_def), parameter :: mintra_power = 0.046822_r_def
    ! Limits of the search for crossings of the curve
    real(kind=r_def), parameter :: min_pressure = 5000.0_r_def
    real(kind=r_def), parameter :: max_height   = 18000.0_r_def

    ! Internal variables
    integer(kind=i_def) :: k, n
    integer(kind=i_def) :: k_limit(2)
    real(kind=r_def)    :: p_limit(2)
    real(kind=r_def)    :: temperature, temperature_below, pressure
    real(kind=r_def)    :: mintra, mintra_below, lapse

    k_limit(:) = 0_i_def
    p_limit(:) = rmdi

    ! Contrails forming at the lowest level: lower limit at the surface
    temperature = theta(map_wth(1) + 1) * exner(map_wth(1) + 1)
    pressure = p_zero * exner(map_wth(1) + 1)**(1.0_r_def / kappa)
    mintra = mintra_coeff * pressure**mintra_power - temperature
    if (mintra > 0.0_r_def) then
      k_limit(1) = 1_i_def
      p_limit(1) = p_zero * exner(map_wth(1))**(1.0_r_def / kappa)
    end if

    ! Find the levels just above the first and last crossings of the curve
    do k = 2, nlayers
      mintra_below = mintra
      temperature = theta(map_wth(1) + k) * exner(map_wth(1) + k)
      pressure = p_zero * exner(map_wth(1) + k)**(1.0_r_def / kappa)
      mintra = mintra_coeff * pressure**mintra_power - temperature
      if ( pressure > min_pressure .and.            &
           height(map_wth(1) + k) < max_height .and. &
           mintra * mintra_below <= 0.0_r_def ) then
        if (k_limit(1) > 0) then
          k_limit(2) = k
        else
          k_limit(1) = k
        end if
      end if
    end do

    ! Pressure of each crossing, assuming a constant lapse rate across the
    ! layer, so that T = T(k) * (p / p(k))**(lapse / g_over_r)
    do n = 1, 2
      k = k_limit(n)
      if (k > 1) then
        temperature = theta(map_wth(1) + k) * exner(map_wth(1) + k)
        temperature_below = theta(map_wth(1) + k - 1) &
                            * exner(map_wth(1) + k - 1)
        pressure = p_zero * exner(map_wth(1) + k)**(1.0_r_def / kappa)
        lapse = (temperature_below - temperature)                    &
                / (height(map_wth(1) + k) - height(map_wth(1) + k - 1))
        p_limit(n) = ( pressure**lapse                                &
                       * (mintra_coeff / temperature)**g_over_r )     &
                     **(1.0_r_def / (lapse - mintra_power * g_over_r))
      end if
    end do

    p_lower(map_2d(1)) = p_limit(1)
    p_upper(map_2d(1)) = p_limit(2)

  end subroutine contrail_limits_code

end module contrail_limits_kernel_mod
