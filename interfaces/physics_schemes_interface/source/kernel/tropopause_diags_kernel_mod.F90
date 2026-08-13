!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Height, temperature, pressure and ICAO height at the tropopause

module tropopause_diags_kernel_mod

use argument_mod,      only : arg_type,          &
                              GH_FIELD, GH_REAL, &
                              GH_READ, GH_WRITE, &
                              CELL_COLUMN,       &
                              ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: tropopause_diags_kernel_type
  private
  ! Args: theta, exner_in_wth, height_wth (in); trop_ht, trop_temp,
  ! trop_press, trop_icao_ht (out)
  type(arg_type) :: meta_args(7) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: tropopause_diags_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: tropopause_diags_code

contains

!> @param[in]     nlayers               Number of layers
!> @param[in]     theta                 Potential temperature
!> @param[in]     exner_in_wth          Exner pressure in wth space
!> @param[in]     height_wth            Height of wth levels above surface
!> @param[in,out] trop_ht               Height at the tropopause
!> @param[in,out] trop_temp             Temperature at the tropopause
!> @param[in,out] trop_press            Pressure at the tropopause
!> @param[in,out] trop_icao_ht          ICAO standard-atmosphere height at
!>                                      the tropopause
!> @param[in]     ndf_wth               No. DOFs per cell for wth space
!> @param[in]     undf_wth              No. unique DOFs for wth space
!> @param[in]     map_wth               Dofmap for wth space column base cell
!> @param[in]     ndf_2d                No. DOFs per cell for 2D space
!> @param[in]     undf_2d               No. unique DOFs for 2D space
!> @param[in]     map_2d                Dofmap for 2D space column base cell
subroutine tropopause_diags_code(nlayers,                    &
                                  theta,                       &
                                  exner_in_wth,                &
                                  height_wth,                  &
                                  trop_ht,                     &
                                  trop_temp,                   &
                                  trop_press,                  &
                                  trop_icao_ht,                &
                                  ndf_wth, undf_wth, map_wth,   &
                                  ndf_2d, undf_2d, map_2d)

  use planet_config_mod, only : p_zero, kappa, gravity, rd

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wth, ndf_2d
  integer(i_def), intent(in) :: undf_wth, undf_2d

  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d

  real(r_def), dimension(undf_wth), intent(in) :: theta, exner_in_wth
  real(r_def), dimension(undf_wth), intent(in) :: height_wth

  real(r_def), dimension(undf_2d), intent(inout) :: trop_ht, trop_temp
  real(r_def), dimension(undf_2d), intent(inout) :: trop_press, trop_icao_ht

  ! Local variables
  integer(i_def) :: k, kk
  integer(i_def) :: lapse_rate_trop_level, cold_point_trop_level
  real(r_def) :: exner_max, exner_min
  real(r_def) :: t_wth(nlayers), lapse_rate(nlayers), lapse_rate_above, dz
  real(r_def) :: g_over_r, lapseupr, lapselwr, delta_lapse
  real(r_def) :: press_at_k, icao_press

  ! Parameters for WMO tropopause definition
  real(r_def), parameter :: lapse_trop = 0.002_r_def   ! K/m
  real(r_def), parameter :: dz_trop = 2000.0_r_def     ! m

  ! Parameters to limit tropopause to given pressure range
  ! (could be set in planet namelist for different planets in future)
  real(r_def), parameter :: p_min_trop = 5000.0_r_def  ! Pa
  real(r_def), parameter :: p_max_trop = 50000.0_r_def ! Pa

  real(r_def), parameter :: vsmall = 1.0e-6_r_def

  ! ICAO standard-atmosphere constants (pressure -> height, Pa -> kft)
  real(r_def), parameter :: ft2m = 0.3048_r_def
  real(r_def), parameter :: mtokft = 1.0_r_def / (ft2m * 1000.0_r_def)
  real(r_def), parameter :: lapse_rate_l = 6.5e-03_r_def
                                                ! For levels below 11,000 gpm
  real(r_def), parameter :: lapse_rate_u = -1.0e-03_r_def
                                                ! For levels above 11,000 gpm
  real(r_def), parameter :: press_bot = 101325.0_r_def
                                                ! ICAO std: surface pressure
  real(r_def), parameter :: press_mid = 22632.0_r_def
                                                !      pressure @ 11,000 gpm
  real(r_def), parameter :: press_top = 5474.87_r_def
                                                !      pressure @ 20,000 gpm
  real(r_def), parameter :: temp_bot = 288.15_r_def ! Surface temperature
  real(r_def), parameter :: temp_top = 216.65_r_def
                                                ! Temperature of isothermal
                                                ! layer
  real(r_def), parameter :: gpm1 = 11000.0_r_def
                                                ! Ht limit (gpm) for std
                                                ! lower lapse rate
  real(r_def), parameter :: gpm2 = 20000.0_r_def
                                                ! Ht (gpm) of top of
                                                ! isothermal layer
  real(r_def) :: zp1, zp2 ! Exponents used for the ICAO calculation

  g_over_r = gravity / rd
  zp1 = lapse_rate_l / g_over_r
  zp2 = lapse_rate_u / g_over_r

  exner_min = (p_min_trop / p_zero)**kappa
  exner_max = (p_max_trop / p_zero)**kappa
  lapse_rate_trop_level = 0
  ! Fallback level must leave room for the k-1 and k+2 accesses used by the
  ! interpolation below; 3 is the lowest level the search loop considers.
  cold_point_trop_level = 3
  t_wth(1) = theta(map_wth(1) + 1) * exner_in_wth(map_wth(1) + 1)
  do k = 2, nlayers
    t_wth(k) = theta(map_wth(1) + k) * exner_in_wth(map_wth(1) + k)
    lapse_rate(k) = (t_wth(k - 1) - t_wth(k))                             &
                  / (height_wth(map_wth(1) + k) -                        &
                     height_wth(map_wth(1) + k - 1))
  end do

  ! Locate the WMO-definition tropopause model level, following the same
  ! pressure-band search as locate_tropopause_kernel_mod. The upper bound
  ! is nlayers - 2 (rather than nlayers - 1) so the k+2 level used by the
  ! interpolation below always stays within the column.
  do k = 3, nlayers - 2
    if (exner_in_wth(map_wth(1) + k - 1) > exner_min .and. &
        exner_in_wth(map_wth(1) + k)   < exner_max) then
      if (t_wth(k) < t_wth(cold_point_trop_level)) then
        ! Set the coldest level to use as a fallback if the lapse-rate
        ! criteria are not met.
        cold_point_trop_level = k
      end if
      if (lapse_rate(k)   < lapse_trop .and. &
          lapse_rate(k - 1) > 0.0_r_def) then
        ! Lapse rate has dropped below the threshold. If this is maintained
        ! for 2km above then the WMO criteria for the tropopause has been
        ! met.
        do kk = k + 1, nlayers
          dz = height_wth(map_wth(1) + kk) - height_wth(map_wth(1) + k)
          if (dz >= dz_trop .or. kk == nlayers) then
            lapse_rate_above = (t_wth(k) - t_wth(kk)) / dz
            exit
          end if
        end do
        if (lapse_rate_above < lapse_trop) then
          lapse_rate_trop_level = k
          exit
        end if
      end if
    end if
  end do

  if (lapse_rate_trop_level > 0) then
    k = lapse_rate_trop_level
  else
    k = cold_point_trop_level
  end if

  ! Lapse rate for the interval below (k-1 ~ k) and above (k+1 ~ k+2) the
  ! tropopause level, used to interpolate height/temperature/pressure at the
  ! crossing point between the two lapse-rate lines.
  lapselwr = lapse_rate(k)
  lapseupr = (t_wth(k + 1) - t_wth(k + 2)) &
           / (height_wth(map_wth(1) + k + 2) - height_wth(map_wth(1) + k + 1))

  delta_lapse = lapselwr - lapseupr
  if (abs(delta_lapse) < vsmall) then
    if (delta_lapse >= 0.0_r_def) delta_lapse = vsmall
    if (delta_lapse <  0.0_r_def) delta_lapse = -vsmall
  end if

  trop_ht(map_2d(1)) =                                                     &
    ((t_wth(k)     + (lapselwr * height_wth(map_wth(1) + k)))   -          &
     (t_wth(k + 1) + (lapseupr * height_wth(map_wth(1) + k + 1)))) / delta_lapse

  if (trop_ht(map_2d(1)) < height_wth(map_wth(1) + k)) then
    ! ensure trop height doesn't undershoot
    trop_ht(map_2d(1)) = height_wth(map_wth(1) + k)
  end if
  if (trop_ht(map_2d(1)) > height_wth(map_wth(1) + k + 1)) then
    ! or overshoot
    trop_ht(map_2d(1)) = height_wth(map_wth(1) + k + 1)
  end if

  trop_temp(map_2d(1)) = t_wth(k) &
    - lapselwr * (trop_ht(map_2d(1)) - height_wth(map_wth(1) + k))

  if (abs(lapselwr) < vsmall) then
    if (lapselwr >= 0.0_r_def) lapselwr = vsmall
    if (lapselwr <  0.0_r_def) lapselwr = -vsmall
  end if

  ! Pressure at the tropopause is derived from the hydrostatic equation.
  press_at_k = p_zero * exner_in_wth(map_wth(1) + k)**(1.0_r_def / kappa)
  trop_press(map_2d(1)) = press_at_k &
    * (trop_temp(map_2d(1)) / t_wth(k))**(g_over_r / lapselwr)

  ! Convert tropopause pressure to ICAO standard-atmosphere height (kft).
  icao_press = trop_press(map_2d(1))
  if (icao_press <= 1000.0_r_def .and. icao_press >= 0.0_r_def) then
    icao_press = 1000.0_r_def
  end if
  if (icao_press > press_bot) then
    icao_press = press_bot
  end if

  if (icao_press > press_mid) then ! Hts up to 11,000 GPM
    icao_press = icao_press / press_bot
    icao_press = 1.0_r_def - icao_press**zp1
    trop_icao_ht(map_2d(1)) = icao_press * temp_bot / lapse_rate_l
  else if (icao_press > press_top) then ! Hts between 11,000 and 20,000 GPM
    icao_press = icao_press / press_mid
    icao_press = -log(icao_press)
    trop_icao_ht(map_2d(1)) = gpm1 + icao_press * temp_top / g_over_r
  else ! Hts above 20,000 GPM
    icao_press = icao_press / press_top
    icao_press = 1.0_r_def - icao_press**zp2
    trop_icao_ht(map_2d(1)) = gpm2 + icao_press * temp_top / lapse_rate_u
  end if
  trop_icao_ht(map_2d(1)) = trop_icao_ht(map_2d(1)) * mtokft

end subroutine tropopause_diags_code

end module tropopause_diags_kernel_mod
