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
                                ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod,    only: Wtheta
use constants_mod,        only: r_def, i_def
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
  type(arg_type) :: meta_args(4) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                   & ! theta
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                   & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                   & ! height_wth
       arg_type(GH_FIELD, GH_INTEGER,GH_WRITE,ANY_DISCONTINUOUS_SPACE_1)& ! trop_level
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
  REAL(KIND=real_umphys) :: lapseupr, lapselwr, delta_lapse

  REAL(KIND=real_umphys), PARAMETER :: vsmall = 1.0e-6
  

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




  ! calculate tropopause height, temperature and pressure, and icao heights, for diagnostics

  ! Assuming the tropopause was found in the steps above, retrieve the
  ! level index for it here...
  ! We shall calculate the lapse rates above and below expected level
  ! and assume linear cross over point is nearer to the true height of the
  ! tropopause level -- it need not correspond to an exact model level --
  ! if no model level was identified return mdi
  
  ! todo: flags for these? (tropopause_height, tropopause_temp, tropopause_pressure)
  ! if level found  
  if (lapse_rate_trop_level > 0) then

    k = lapse_rate_trop_level

    ! lapse rate for interval above, k+1 -> k+2
    ! lapse rate for interval below, k-1 -> k
    lapse_upr = lapse_rate(k+2)
    lapse_lwr = lapse_rate(k)

    delta_lapse = lapselwr - lapseupr
    IF ( ABS(delta_lapse) < vsmall ) THEN
      IF ( delta_lapse >= 0 ) delta_lapse =  vsmall
      IF ( delta_lapse <  0 ) delta_lapse = -vsmall
    END IF

    ! height of tropopause
    ! note: the um code subtracts planet_radius from height
    tropopause_height(map_2d(1)) = (                       &
      (t_wth(k) + (lapselwr*height_wth(k)))     &
    - (t_wth(k+1) + (lapseupr*height_wth(k+1))) &
    ) / delta_lapse

    ! temperature at tropopause
    tropopause_temp(map_2d(1)) = t_wth(k) -                              &
      lapselwr * (tropopause_height(i,j) - height_wth(i,j,k))

    ! pressure at tropopause is derived from hydrostatic equation
    ! todo: we'll need an equivalent of p_theta_levels
    IF ( ABS(lapselwr) < vsmall ) THEN
      IF ( lapseLwr >= 0 ) lapselwr =  vsmall
      IF ( lapselwr <  0 ) lapselwr = -vsmall
    END IF
    
    tropopause_press(map_2d(1)) = p_theta_levels(k) *                       &
      (tropopause_temp(map_2d(1)) / temp(k)) ** (g_over_r/lapselwr)

  endif


  ! calculate tropopauise icao heights, temperature and pressure



end subroutine locate_tropopause_code

end module locate_tropopause_kernel_mod
