!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Tropopause aviation diagnostics

module trop_diags_kernel_mod

use argument_mod,      only : arg_type,                       &
                              GH_FIELD, GH_REAL, GH_SCALAR,   &
                              GH_READ, GH_WRITE, CELL_COLUMN, &
                              ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

type, public, extends(kernel_type) :: trop_diags_kernel_type
  private
  type(arg_type) :: meta_args(12) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! theta
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! height_wth
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! g_over_r
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_press_flag
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_press
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_temp_flag
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_temp
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_height_flag
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_height
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_icao_height_flag
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! trop_icao_height
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: trop_diags_code
end type

public :: trop_diags_code

contains

!> @brief Calculate tropopause aviation diagnostics (pressure, temperature and height)
!> @details If a field is not empty, its dependencies are guaranteed by the caller.
!>          If the tropopause is not found, outputs are set to missing data.
!> @param[in]  nlayers                Number of layers
!> @param[in]  theta                  Potential temperature
!> @param[in]  exner_in_wth           Exner pressure in wth space
!> @param[in]  height_wth             Height of wth levels above surface
!> @param[in]  g_over_r               Gravity / specific dry air gas constant
!> @param[out] trop_press             Pressure at the tropopause
!> @param[out] trop_temp              Temperature at the tropopause
!> @param[out] trop_height            Height of the tropopause above surface
!> @param[out] trop_icao_height       ICAO height of the tropopause
!> @param[in]  ndf_wth                No. DOFs per cell for wth space
!> @param[in]  undf_wth               No. unique DOFs for wth space
!> @param[in]  map_wth                Dofmap for wth space column base cell
!> @param[in]  ndf_2d                 No. DOFs per cell for 2D space
!> @param[in]  undf_2d                No. unique DOFs for 2D space
!> @param[in]  map_2d                 Dofmap for 2D space column base cell
subroutine trop_diags_code(nlayers,                                  &
                           theta,                                    &
                           exner_in_wth,                             &
                           height_wth,                               &
                           g_over_r,                                 &
                           trop_pres_flag, trop_press,               &
                           trop_temp_flag, trop_temp,                &
                           trop_height_flag, trop_height,            &
                           trop_icao_height_flag, trop_icao_height,  &
                           ndf_wth, undf_wth, map_wth,               &
                           ndf_2d, undf_2d, map_2d)

  use planet_config_mod,        only : p_zero, kappa
  use missing_data_mod,         only : rmdi
  use icao_heights_kernel_mod,  only : icao_heights_kernel_code
  use empty_data_mod,           only : empty_real_data

  implicit none

  ! Arguments (algorithm)
  real(r_def), dimension(undf_wth), intent(in) :: theta
  real(r_def), dimension(undf_wth), intent(in) :: exner_in_wth
  real(r_def), dimension(undf_wth), intent(in) :: height_wth
  real(r_def), intent(in)                      :: g_over_r

  logical(l_def), intent(in)                   :: trop_pres_flag
  real(r_def), dimension(undf_2d), intent(out) :: trop_press(:)

  logical(l_def), intent(in)                   :: trop_temp_flag
  real(r_def), dimension(undf_2d), intent(out) :: trop_temp(:)

  logical(l_def), intent(in)                   :: trop_height_flag
  real(r_def), dimension(undf_2d), intent(out) :: trop_height(:)

  logical(l_def), intent(in)                   :: trop_icao_height_flag
  real(r_def), dimension(undf_2d), intent(out) :: trop_icao_height(:)

  ! Arguments (kernel)
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wth, ndf_2d
  integer(i_def), intent(in) :: undf_wth, undf_2d

  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d


  ! Local variables
  integer(i_def) :: k, kk
  integer(i_def) :: lapse_rate_trop_level
  real(r_def) :: exner_max, exner_min
  real(r_def) :: t_wth(nlayers), lapse_rate(nlayers), lapse_rate_above, dz
  real(r_def) :: lapse_upr, lapse_lwr, delta_lapse, press_wth

  ! Parameters for WMO tropopause definition
  real(r_def), parameter :: lapse_trop = 0.002_r_def   ! K/m
  real(r_def), parameter :: dz_trop = 2000.0_r_def     ! m

  ! Parameters to limitguaranteed tropopause to given pressure range
  real(r_def), parameter :: p_min_trop = 5000.0_r_def  ! Pa
  real(r_def), parameter :: p_max_trop = 50000.0_r_def ! Pa

  ! Small number used to protect divisions from near-zero denominators
  real(r_def), parameter :: vsmall = 1.0e-9_r_def


  exner_min = (p_min_trop/p_zero)**kappa
  exner_max = (p_max_trop/p_zero)**kappa

  ! Initialise all outputs to missing data. They are only meaningful when a
  ! lapse-rate (WMO) tropopause has been located.
  if (trop_pres_flag)        trop_press(map_2d(1))       = rmdi
  if (trop_temp_flag)        trop_temp(map_2d(1))        = rmdi
  if (trop_height_flag)      trop_height(map_2d(1))      = rmdi
  if (trop_icao_height_flag) trop_icao_height(map_2d(1)) = rmdi

  ! Temperature and lapse rate on wth levels
  lapse_rate_trop_level = 0
  t_wth(1) = theta(map_wth(1)+1) * exner_in_wth(map_wth(1)+1)
  do k=2, nlayers
    t_wth(k) = theta(map_wth(1)+k) * exner_in_wth(map_wth(1)+k)
    lapse_rate(k) = ( t_wth(k-1) - t_wth(k) ) &
                  / ( height_wth(map_wth(1)+k) - height_wth(map_wth(1)+k-1) )
  end do

  ! Locate the lapse-rate (WMO) tropopause
  do k=3, nlayers-1
    if (exner_in_wth(map_wth(1)+k-1) > exner_min .and. &
        exner_in_wth(map_wth(1)+k)   < exner_max) then
      if (lapse_rate(k)   < lapse_trop .and. &
          lapse_rate(k-1) > 0.0_r_def) then
        ! Lapse rate has dropped below the threshold. If this is maintained
        ! for 2km above then the WMO criteria for the tropopause has been met.
        do kk=k+1, nlayers
          dz = height_wth(map_wth(1)+kk) - height_wth(map_wth(1)+k)
          if (dz >= dz_trop .or. kk==nlayers) then
            lapse_rate_above = ( t_wth(k) - t_wth(kk) ) / dz
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


  ! if level found
  if (lapse_rate_trop_level > 0) then

    k = lapse_rate_trop_level

    ! lapse rate for interval above, k+1 -> k+2
    ! lapse rate for interval below, k-1 -> k
    lapse_upr = lapse_rate(min(k+2, nlayers))
    lapse_lwr = lapse_rate(k)

    delta_lapse = lapse_lwr - lapse_upr
    if ( abs(delta_lapse) < vsmall ) then
      if ( delta_lapse >= 0.0_r_def ) delta_lapse =  vsmall
      if ( delta_lapse <  0.0_r_def ) delta_lapse = -vsmall
    end if

    ! height of tropopause between k and k+1
    if (trop_height_flag) then
      trop_height(map_2d(1)) = (                                     &
          (t_wth(k)   + (lapse_lwr * height_wth(map_wth(1)+k)))      &
        - (t_wth(k+1) + (lapse_upr * height_wth(map_wth(1)+k+1)))    &
        ) / delta_lapse
    end if

    ! temperature at tropopause
    if (trop_temp_flag) then
      trop_temp(map_2d(1)) = t_wth(k) -                              &
        lapse_lwr * (trop_height(map_2d(1)) - height_wth(map_wth(1)+k))
    end if

    ! pressure at tropopause is derived from the hydrostatic equation
    if (trop_pres_flag) then
      if ( abs(lapse_lwr) < vsmall ) then
        if ( lapse_lwr >= 0.0_r_def ) lapse_lwr =  vsmall
        if ( lapse_lwr <  0.0_r_def ) lapse_lwr = -vsmall
      end if
      press_wth = p_zero * exner_in_wth(map_wth(1)+k) ** (1.0_r_def / kappa)
      trop_press(map_2d(1)) = press_wth *                           &
        (trop_temp(map_2d(1)) / t_wth(k)) ** (g_over_r / lapse_lwr)
    end if

    ! ICAO height of the tropopause
    if (trop_icao_height_flag) then
      call icao_heights_kernel_code( &
        nlayers, trop_icao_height, trop_press, g_over_r, &
        ndf_2d, undf_2d, map_2d)
    end if

  end if

end subroutine trop_diags_code

end module trop_diags_kernel_mod

