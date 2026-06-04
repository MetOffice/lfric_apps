!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Locate tropopause level

module locate_tropopause_kernel_mod

use argument_mod,         only: arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_READ, GH_WRITE, &
                                CELL_COLUMN,       &
                                GH_INTEGER,        &
                                GH_SCALAR,         &
                                GH_LOGICAL,        &
                                ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod,    only: Wtheta
use constants_mod,        only: r_def, i_def, rmdi
use kernel_mod,           only: kernel_type
use planet_constants_mod, only: g_over_r

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: locate_tropopause_kernel_type
  private
  type(arg_type) :: meta_args(12) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! theta
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! height_wth
       arg_type(GH_FIELD, GH_INTEGER,GH_WRITE,ANY_DISCONTINUOUS_SPACE_1),& ! trop_level
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_height_flag
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_temp_flag
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_press_flag
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_icao_height_flag
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_height
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_temp
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_press
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! trop_icao_height
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: locate_tropopause_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: locate_tropopause_code

contains

!> @param[in]     nlayers               Number of layers
!> @param[in]     theta                 Potential temperature
!> @param[in]     exner_in_wth          Exner pressure in wth space
!> @param[in]     height_wth            Height of wth levels above surface
!> @param[in,out] trop_level            Level of tropopause
!> @param[in]     trop_height_flag      Whether to calculate tropopause height
!> @param[in]     trop_temp_flag        Whether to calculate tropopause temperature
!> @param[in]     trop_press_flag       Whether to calculate tropopause pressure
!> @param[in]     trop_icao_height_flag Whether to calculate tropopause icao height
!> @param[in,out] trop_height           Height of tropopause
!> @param[in,out] trop_temp             Temperature at tropopause
!> @param[in,out] trop_press            Pressure at tropopause
!> @param[in,out] trop_icao_height      ICAO height of tropopause
!> @param[in]     ndf_wth               No. DOFs per cell for wth space
!> @param[in]     undf_wth              No. unique DOFs for wth space
!> @param[in]     map_wth               Dofmap for wth space column base cell
!> @param[in]     ndf_2d                No. DOFs per cell for 2D space
!> @param[in]     undf_2d               No. unique DOFs for 2D space
!> @param[in]     map_2d                Dofmap for 2D space column base cell
subroutine locate_tropopause_code(nlayers,                    &
                                  theta,                      &
                                  exner_in_wth,               &
                                  height_wth,                 &
                                  trop_level,                 &
                                  trop_height_flag,           &
                                  trop_temp_flag,             &
                                  trop_press_flag,            &
                                  trop_icao_height_flag,      &
                                  trop_height,                &
                                  trop_temp,                  &
                                  trop_press,                 &
                                  trop_icao_height,           &
                                  ndf_wth, undf_wth, map_wth, &
                                  ndf_2d, undf_2d, map_2d)

  use planet_config_mod, only : p_zero, kappa

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wth, ndf_2d
  integer(i_def), intent(in) :: undf_wth, undf_2d

  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d

  real(r_def), dimension(undf_wth), intent(in) :: theta, exner_in_wth
  real(r_def), dimension(undf_wth), intent(in) :: height_wth

  integer(i_def), dimension(undf_2d), intent(inout) :: trop_level
  real(r_def), dimension(undf_2d), intent(inout) :: trop_height
  real(r_def), dimension(undf_2d), intent(inout) :: trop_temp
  real(r_def), dimension(undf_2d), intent(inout) :: trop_press
  real(r_def), dimension(undf_2d), intent(inout) :: trop_icao_height

  logical, intent(in) :: trop_height_flag, trop_temp_flag, trop_press_flag, trop_icao_height_flag

  ! Local variables
  integer(i_def) :: k, kk
  integer(i_def) :: lapse_rate_trop_level, cold_point_trop_level
  real(r_def) :: exner_max, exner_min
  real(r_def) :: t_wth(nlayers), lapse_rate(nlayers), lapse_rate_above, dz

  ! Parameters for WMO tropopause definition
  real(r_def), parameter :: lapse_trop = 0.002_r_def   ! K/m
  real(r_def), parameter :: dz_trop = 2000.0_r_def     ! m

  ! Parameters to limit tropopause to given pressure range
  ! (could be set in planet namelist for different planets in future)
  real(r_def), parameter :: p_min_trop = 5000.0_r_def  ! Pa
  real(r_def), parameter :: p_max_trop = 50000.0_r_def ! Pa

  ! Lapse rates below and above current layer, and the difference between them
  REAL(r_def) :: lapse_upr, lapse_lwr, delta_lapse

  REAL(r_def), PARAMETER :: vsmall = 1.0e-6
  REAL(r_def) :: press_wth
  

  exner_min = (p_min_trop/p_zero)**kappa
  exner_max = (p_max_trop/p_zero)**kappa
  lapse_rate_trop_level = 0
  cold_point_trop_level = 1
  t_wth(1) = theta(map_wth(1)+1) * exner_in_wth(map_wth(1)+1)
  do k=2, nlayers
    t_wth(k) = theta(map_wth(1)+k) * exner_in_wth(map_wth(1)+k)
    lapse_rate(k) = ( t_wth(k-1) - t_wth(k) ) &
                  / ( height_wth(map_wth(1)+k) - height_wth(map_wth(1)+k-1) )
  end do
  do k=3, nlayers-1
    if (exner_in_wth(map_wth(1)+k-1) > exner_min .and. &
        exner_in_wth(map_wth(1)+k)   < exner_max) then
      if (t_wth(k) < t_wth(cold_point_trop_level)) then
        ! Set the coldest level to use as a fallback if the lapse-rate
        ! criteria are not met.
        cold_point_trop_level = k
      end if
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
  if (lapse_rate_trop_level > 0) then
    trop_level(map_2d(1)) = lapse_rate_trop_level
  else
    trop_level(map_2d(1)) = cold_point_trop_level
  end if




  ! Calculate tropopause height, temperature and pressure, and icao heights
  ! for diagnostics, from the tropopause level identified above.
  !
  ! We shall calculate the lapse rates above and below expected level
  ! and assume linear cross over point is nearer to the true height of the
  ! tropopause level -- it need not correspond to an exact model level --
  ! if no model level was identified return mdi
  
  ! if level found  
  if (lapse_rate_trop_level > 0) then

    k = lapse_rate_trop_level

    ! lapse rate for interval above, k+1 -> k+2
    ! lapse rate for interval below, k-1 -> k
    lapse_upr = lapse_rate(k+2)
    lapse_lwr = lapse_rate(k)

    delta_lapse = lapse_lwr - lapse_upr
    IF ( ABS(delta_lapse) < vsmall ) THEN
      IF ( delta_lapse >= 0 ) delta_lapse =  vsmall
      IF ( delta_lapse <  0 ) delta_lapse = -vsmall
    END IF

    ! height of tropopause between k and k+1
    if (trop_height_flag) then
      ! discuss: the um code subtracts planet_radius from height
      trop_height(map_2d(1)) = (                         &
        (t_wth(k) + (lapse_lwr * height_wth(map_wth(1)+k)))     &
      - (t_wth(k+1) + (lapse_upr * height_wth(map_wth(1)+k+1))) &
      ) / delta_lapse
    end if

    ! temperature at tropopause
    if (trop_temp_flag) then
      trop_temp(map_2d(1)) = t_wth(k) -                              &
        lapse_lwr * (trop_height(map_2d(1)) - height_wth(map_wth(1) + k))
    end if

    ! pressure at tropopause is derived from hydrostatic equation
    if (trop_press_flag) then
      IF ( ABS(lapse_lwr) < vsmall ) THEN
        IF ( lapse_lwr >= 0 ) lapse_lwr =  vsmall
        IF ( lapse_lwr <  0 ) lapse_lwr = -vsmall
      END IF   

      press_wth = p_zero * exner_in_wth(map_wth(1)+k) ** (1.0_r_def / kappa)
      trop_press(map_2d(1)) = press_wth *                             &
        (trop_temp(map_2d(1)) / t_wth(k)) ** (g_over_r/lapse_lwr)
    end if

  endif


  ! calculate tropopause icao height of the tropopause 
  ! todo: call the kernel from https://github.com/MetOffice/lfric_apps/pull/169
  if (trop_icao_height_flag) then
    ! call diags_icao_heights_kernel_code( &
    !   nlayers, trop_icao_height, trop_press, ndf_2d, undf_2d, map_2d)
    trop_icao_height(map_2d(1)) = rmdi
  end if



end subroutine locate_tropopause_code

end module locate_tropopause_kernel_mod
