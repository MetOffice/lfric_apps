!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Tropopause aviation diagnostics

module trop_diags_kernel_mod

use argument_mod,      only : arg_type,                                  &
                              GH_FIELD, GH_REAL, GH_SCALAR, GH_LOGICAL,  &
                              GH_READ, GH_WRITE, CELL_COLUMN,            &
                              ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

! layered input fields before flat output fields for nlayers > 1
type, public, extends(kernel_type) :: trop_diags_kernel_type
  private
  type(arg_type) :: meta_args(12) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! theta_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! exner_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! height_wth
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! g_over_r
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_press
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_temp
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_height
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! trop_icao_height
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_press_flag
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_temp_flag
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                         & ! trop_height_flag
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                          & ! trop_icao_height_flag
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
!> @param[in]  theta_wth              Potential temperature
!> @param[in]  exner_wth              Exner pressure in wth space
!> @param[in]  height_wth             Height of wth levels above surface
!> @param[in]  g_over_r               Gravity / specific dry air gas constant
!> @param[out] trop_press             Pressure at the tropopause
!> @param[out] trop_temp              Temperature at the tropopause
!> @param[out] trop_height            Height of the tropopause above surface
!> @param[out] trop_icao_height       ICAO height of the tropopause
!> @param[in]  trop_press_flag        Request flag
!> @param[in]  trop_temp_flag         Request flag
!> @param[in]  trop_height_flag       Request flag
!> @param[in]  trop_icao_height_flag  Request flag
!> @param[in]  ndf_wth                No. DOFs per cell for wth space
!> @param[in]  undf_wth               No. unique DOFs for wth space
!> @param[in]  map_wth                Dofmap for wth space column base cell
!> @param[in]  ndf_2d                 No. DOFs per cell for 2D space
!> @param[in]  undf_2d                No. unique DOFs for 2D space
!> @param[in]  map_2d                 Dofmap for 2D space column base cell
subroutine trop_diags_code(nlayers,                    &
                           theta_wth,                  &
                           exner_wth,                  &
                           height_wth,                 &
                           g_over_r,                   &
                           trop_press,                 &
                           trop_temp,                  &
                           trop_height,                &
                           trop_icao_height,           &
                           trop_pres_flag,             &
                           trop_temp_flag,             &
                           trop_height_flag,           &
                           trop_icao_height_flag,      &
                           ndf_wth, undf_wth, map_wth, &
                           ndf_2d, undf_2d, map_2d)

  use planet_config_mod,        only : p_zero, kappa
  use missing_data_mod,         only : rmdi, imdi
  use icao_heights_kernel_mod,  only : icao_heights_kernel_code
  use empty_data_mod,           only : empty_real_data

  use log_mod,                        only: log_event,       &
                                            LOG_LEVEL_ALWAYS,  &
                                            log_scratch_space

  implicit none

  ! Arguments (algorithm)
  real(r_def), dimension(undf_wth), intent(in) :: theta_wth
  real(r_def), dimension(undf_wth), intent(in) :: exner_wth
  real(r_def), dimension(undf_wth), intent(in) :: height_wth
  real(r_def), intent(in)                      :: g_over_r

  real(r_def), dimension(undf_2d), intent(out) :: trop_press(:)
  real(r_def), dimension(undf_2d), intent(out) :: trop_temp(:)
  real(r_def), dimension(undf_2d), intent(out) :: trop_height(:)
  real(r_def), dimension(undf_2d), intent(out) :: trop_icao_height(:)

  logical(l_def), intent(in)                   :: trop_pres_flag
  logical(l_def), intent(in)                   :: trop_temp_flag
  logical(l_def), intent(in)                   :: trop_height_flag
  logical(l_def), intent(in)                   :: trop_icao_height_flag


  ! Arguments (kernel)
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_2d, ndf_wth
  integer(i_def), intent(in) :: undf_2d, undf_wth

  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d
  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth


  ! Local variables
  integer(i_def) :: k, k2km
  integer(i_def) :: trop_level
  real(r_def) :: t_wth(nlayers)
  real(r_def) :: lapse, lapse_below, lapse_above, lapse_2km
  real(r_def) :: dz
  real(r_def) :: delta_lapse, press_wth

  ! Parameters for WMO tropopause definition
  real(r_def), parameter :: lapse_trop = 0.002_r_def   ! K/m
  real(r_def), parameter :: dz_trop = 2000.0_r_def     ! m

  ! Parameters to limit tropopause to given pressure range
  real(r_def), parameter :: p_min_trop = 5000.0_r_def  ! Pa
  real(r_def), parameter :: p_max_trop = 50000.0_r_def ! Pa

  ! Small number used to protect divisions from near-zero denominators
  real(r_def), parameter :: vsmall = 1.0e-9_r_def


  ! cut off limits to be used in tropopause calculations.
  ! todo: put this in a constants module? in the um it was in pws_diags_mod
  ! arbritary limits for high and low trop levels for search
  real(r_def), parameter :: heightcut_top = 22000.0
  real(r_def), parameter :: heightcut_bot = 4500.0
  ! max temp allowed for tropopause
  real(r_def), parameter :: tempcut = 243.0


  logical(l_def), save :: print_once = .true.


  ! Initialise all outputs to missing data. They are only meaningful when a
  ! lapse-rate (WMO) tropopause has been located.
  ! todo: do this to fields which are initialised (not just the requested)
  if (trop_pres_flag)        trop_press(map_2d(1))       = rmdi
  if (trop_temp_flag)        trop_temp(map_2d(1))        = rmdi
  if (trop_height_flag)      trop_height(map_2d(1))      = rmdi
  if (trop_icao_height_flag) trop_icao_height(map_2d(1)) = rmdi

  ! Temperature on wth levels
  ! todo: put into following loop
  do k=1, nlayers
    t_wth(k) = theta_wth(map_wth(1)+k) * exner_wth(map_wth(1)+k)
  end do

  if (print_once) then
     write(log_scratch_space, *) 'nlayers=', nlayers
     call log_event(log_scratch_space, LOG_LEVEL_ALWAYS)
  end if


  ! Locate the lapse-rate (WMO) tropopause
  trop_level = imdi
  do k=1, nlayers

    if (t_wth(k) < tempcut .and.             &
        height_wth(k) > heightcut_bot .and.  &
        height_wth(k) < heightcut_top) then

      lapse = (t_wth(k) - t_wth(k+1)) / (height_wth(k+1) - height_wth(k))
      lapse_below = (t_wth(k-1) - t_wth(k)) / (height_wth(k) - height_wth(k-1))

      if (lapse < lapse_trop .and. lapse_below > 0.0_r_def) then
        ! Lapse rate has dropped below the threshold. If this is maintained
        ! for 2km above then the WMO criteria for the tropopause has been met.
        do k2km=k, nlayers
          if (height_wth(k2km) > heightcut_top) exit

          if ((height_wth(k2km)-height_wth(k)) >= 2000.0 ) then
            lapse_2km = (t_wth(k) - t_wth(k2km)) / (height_wth(k2km) - height_wth(k))
            ! if 2km interval also < 2 then we have the tropopause level
            if (lapse_2km < lapse_trop) then
              trop_level = k
              exit
            end if
          end if

        end do  ! looking upwards for 2km

      end if  ! lapse values in range
    end if
  end do

  if (print_once) then
     write(log_scratch_space, '(A, I0)') 'trop_level=', trop_level
     call log_event(log_scratch_space, LOG_LEVEL_ALWAYS)
  end if


  ! if level found
  if (trop_level > 0) then
    k = trop_level

    ! todo: lapse_below was already calculated in the loop above
    lapse_above = (t_wth(k+1) - t_wth(k+2)) / (height_wth(k+2) - height_wth(k+1))
    lapse_below = (t_wth(k-1) - t_wth(k)) / (height_wth(k) - height_wth(k-1))
    delta_lapse = lapse_below - lapse_above
    if ( abs(delta_lapse) < vsmall ) then
      if ( delta_lapse >= 0.0_r_def ) delta_lapse =  vsmall
      if ( delta_lapse <  0.0_r_def ) delta_lapse = -vsmall
    end if

    ! height of tropopause between k and k+1
    if (trop_height_flag) then
      trop_height(map_2d(1)) = (                                     &
          (t_wth(k)   + (lapse_below * height_wth(map_wth(1)+k)))  &
        - (t_wth(k+1) + (lapse_above * height_wth(map_wth(1)+k+1)))    &
        ) / delta_lapse

      ! ensure trop level doesn't undershoot
      if (trop_height(map_2d(1)) < height_wth(map_wth(1)+k)) then
        trop_height(map_2d(1)) = height_wth(map_wth(1)+k)
      end if
      ! or overshoot
      if (trop_height(map_2d(1)) > height_wth(map_wth(1)+k+1)) then
        trop_height(map_2d(1)) = height_wth(map_wth(1)+k+1)
      end if

    end if

    ! temperature at tropopause
    if (trop_temp_flag) then
      trop_temp(map_2d(1)) = t_wth(k) -                                     &
        lapse_below * (trop_height(map_2d(1)) - height_wth(map_wth(1)+k-1))
    end if

    ! pressure at tropopause is derived from the hydrostatic equation
    if ( abs(lapse_below) < vsmall ) then
      if ( lapse_below >= 0.0_r_def ) lapse_below =  vsmall
      if ( lapse_below <  0.0_r_def ) lapse_below = -vsmall
    end if
    if (trop_pres_flag) then
      press_wth = p_zero * exner_wth(map_wth(1)+k) ** (1.0_r_def / kappa)
      trop_press(map_2d(1)) = press_wth *                           &
        (trop_temp(map_2d(1)) / t_wth(k)) ** (g_over_r / lapse_below)
    end if

    ! ICAO height of the tropopause
    if (trop_icao_height_flag) then
      call icao_heights_kernel_code( &
        nlayers, trop_icao_height, trop_press, g_over_r, &
        ndf_2d, undf_2d, map_2d)
    end if

  end if

  print_once = .false.

end subroutine trop_diags_code

end module trop_diags_kernel_mod

