!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Relative humidity on model levels for the icing potential
!>        diagnostic.

module icing_pot_rh_kernel_mod

  use argument_mod,  only: arg_type,            &
                           GH_FIELD, GH_SCALAR, &
                           GH_READ, GH_WRITE,   &
                           GH_REAL, CELL_COLUMN
  use constants_mod, only: r_def, i_def
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: icing_pot_rh_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta), & ! rh
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), & ! temperature
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), & ! exner
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), & ! qv
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! p_zero
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! kappa
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          & ! repsilon
         arg_type(GH_SCALAR, GH_REAL, GH_READ)           & ! zerodegc
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: icing_pot_rh_code
  end type icing_pot_rh_kernel_type

  public :: icing_pot_rh_code

contains

  !> @brief   Calculates relative humidity on model levels for the icing
  !>          potential diagnostic.
  !> @details The saturation vapour pressure uses the Magnus-type
  !>          expression es = 6.11 * 10**(7.5 T / (237.7 + T)), with T in
  !>          degrees Celsius and es in hPa. The saturation mixing ratio is
  !>          then epsilon * es / (p - es), with p in hPa, and the relative
  !>          humidity is the vapour mixing ratio q / (1 - q) divided by the
  !>          saturation mixing ratio.
  !> @param[in]     nlayers     Number of layers
  !> @param[in,out] rh          Relative humidity (fraction)
  !> @param[in]     temperature Temperature (K)
  !> @param[in]     exner       Exner pressure
  !> @param[in]     qv          Specific humidity (kg kg-1)
  !> @param[in]     p_zero      Reference pressure (Pa)
  !> @param[in]     kappa       Rd / cp
  !> @param[in]     repsilon    Ratio of gas constants, Rd / Rv
  !> @param[in]     zerodegc    Zero degrees Celsius in Kelvin
  !> @param[in]     ndf_wth     Number of degrees of freedom per cell
  !> @param[in]     undf_wth    Number of total degrees of freedom
  !> @param[in]     map_wth     Dofmap for the cell at the base of the column
  subroutine icing_pot_rh_code(nlayers,                    &
                               rh,                         &
                               temperature,                &
                               exner,                      &
                               qv,                         &
                               p_zero,                     &
                               kappa,                      &
                               repsilon,                   &
                               zerodegc,                   &
                               ndf_wth, undf_wth, map_wth)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(inout), dimension(undf_wth) :: rh
    real(kind=r_def), intent(in),    dimension(undf_wth) :: temperature
    real(kind=r_def), intent(in),    dimension(undf_wth) :: exner
    real(kind=r_def), intent(in),    dimension(undf_wth) :: qv
    real(kind=r_def), intent(in) :: p_zero, kappa, repsilon, zerodegc

    ! Saturation vapour pressure at zero degrees Celsius (hPa)
    real(kind=r_def), parameter :: es_zero = 6.11_r_def
    ! Multiplier of temperature in the exponent of the Magnus expression
    real(kind=r_def), parameter :: magnus_a = 7.5_r_def
    ! 237.7 - 273.15: Magnus denominator constant (237.7 degrees Celsius)
    ! written relative to temperature in Kelvin
    real(kind=r_def), parameter :: magnus_b = 35.45_r_def
    ! Number of Pa in one hPa
    real(kind=r_def), parameter :: pa_per_hpa = 100.0_r_def
    real(kind=r_def), parameter :: ten = 10.0_r_def

    ! Internal variables
    integer(kind=i_def) :: k
    real(kind=r_def)    :: t, pressure_hpa, es, sat_mr

    do k = 0, nlayers
      t = temperature(map_wth(1) + k)
      pressure_hpa = p_zero * exner(map_wth(1) + k)**(1.0_r_def / kappa) &
                     / pa_per_hpa
      es = es_zero * ten**( magnus_a * (t - zerodegc) / (t - magnus_b) )
      sat_mr = repsilon * es / (pressure_hpa - es)
      rh(map_wth(1) + k) = ( qv(map_wth(1) + k)                          &
                             / (1.0_r_def - qv(map_wth(1) + k)) ) / sat_mr
    end do

  end subroutine icing_pot_rh_code

end module icing_pot_rh_kernel_mod
