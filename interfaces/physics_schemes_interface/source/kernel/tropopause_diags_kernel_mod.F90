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
                              GH_SCALAR,         &
                              CELL_COLUMN,       &
                              ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def, rmdi
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
  ! trop_press, trop_icao_ht (out); p_zero, kappa, g_over_r (in)
  type(arg_type) :: meta_args(10) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR, GH_REAL, GH_READ)                             &
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

!> @brief   Compute tropopause height, temperature, pressure and ICAO height.
!> @details The caller guarantees that whenever trop_press is not empty,
!>          trop_temp is also not empty, and whenever trop_icao_ht is not
!>          empty, trop_press is also not empty.
!> @param[in]     nlayers               Number of layers
!> @param[in]     theta                 Potential temperature
!> @param[in]     exner_in_wth          Exner pressure in wth space
!> @param[in]     height_wth            Height of wth levels above surface
!> @param[in,out] trop_ht               Height at the tropopause
!> @param[in,out] trop_temp             Temperature at the tropopause
!> @param[in,out] trop_press            Pressure at the tropopause
!> @param[in,out] trop_icao_ht          ICAO standard-atmosphere height at
!>                                      the tropopause
!> @param[in]     p_zero                Reference surface pressure
!> @param[in]     kappa                 R/cp for dry air
!> @param[in]     g_over_r              Gravity divided by the dry air gas
!>                                      constant, for the pressure and ICAO
!>                                      height calculations
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
                                  p_zero, kappa, g_over_r,     &
                                  ndf_wth, undf_wth, map_wth,   &
                                  ndf_2d, undf_2d, map_2d)

  use empty_data_mod,          only : empty_real_data
  use icao_heights_kernel_mod, only : icao_heights_kernel_code

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wth, ndf_2d
  integer(i_def), intent(in) :: undf_wth, undf_2d

  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d

  real(r_def), dimension(undf_wth), intent(in) :: theta, exner_in_wth
  real(r_def), dimension(undf_wth), intent(in) :: height_wth

  real(r_def), dimension(undf_2d), intent(inout) :: trop_ht

  ! trop_temp/trop_press/trop_icao_ht may be associated with a shared
  ! empty placeholder (empty_real_data) when not requested; each block
  ! below is skipped in that case rather than writing into it.
  real(r_def), pointer, dimension(:), intent(inout) :: trop_temp
  real(r_def), pointer, dimension(:), intent(inout) :: trop_press
  real(r_def), pointer, dimension(:), intent(inout) :: trop_icao_ht

  real(r_def), intent(in) :: p_zero, kappa, g_over_r

  ! Local variables
  integer(i_def) :: k, kk
  integer(i_def) :: lapse_rate_trop_level
  real(r_def) :: t_wth(nlayers), lapse_rate(nlayers), lapse_rate_above, dz
  real(r_def) :: lapseupr, lapselwr, delta_lapse
  real(r_def) :: press_at_k

  ! Parameters for WMO tropopause definition
  real(r_def), parameter :: lapse_trop = 0.002_r_def   ! K/m
  real(r_def), parameter :: dz_trop = 2000.0_r_def     ! m

  ! Height/temperature band used to gate candidate tropopause levels,
  ! ported verbatim from UM's pws_diags_mod.F90 (heightcut_bot,
  ! heightcut_top, tempcut) and pws_tropoht_mod.F90 (the search using them),
  ! rather than reusing locate_tropopause_kernel_mod's Exner-band criterion,
  ! to keep this diagnostic's output as close to the UM reference as
  ! reasonably possible.
  ! (could be set in planet namelist for different planets in future)
  real(r_def), parameter :: heightcut_bot = 4500.0_r_def  ! m
  real(r_def), parameter :: heightcut_top = 32000.0_r_def ! m
  real(r_def), parameter :: tempcut = 243.0_r_def         ! K

  real(r_def), parameter :: vsmall = 1.0e-6_r_def

  lapse_rate_trop_level = 0
  t_wth(1) = theta(map_wth(1) + 1) * exner_in_wth(map_wth(1) + 1)
  do k = 2, nlayers
    t_wth(k) = theta(map_wth(1) + k) * exner_in_wth(map_wth(1) + k)
    lapse_rate(k) = (t_wth(k - 1) - t_wth(k))                             &
                  / (height_wth(map_wth(1) + k) -                        &
                     height_wth(map_wth(1) + k - 1))
  end do

  ! Locate the WMO-definition tropopause model level, following UM's
  ! height + temperature band search (pws_tropoht_mod.F90) rather than
  ! locate_tropopause_kernel_mod's Exner-band criterion. The upper bound
  ! is nlayers - 2 (rather than nlayers - 1) so the k+2 level used by the
  ! interpolation below always stays within the column.
  do k = 3, nlayers - 2
    if (height_wth(map_wth(1) + k) > heightcut_bot .and. &
        height_wth(map_wth(1) + k) < heightcut_top .and. &
        t_wth(k) < tempcut) then
      ! lapse_rate(k + 1) is the (k ~ k+1) interval and lapse_rate(k) is the
      ! (k-1 ~ k) interval below it - UM's "lapse"/"lapse_below" pair tested
      ! at level k, not lapse_rate(k)/lapse_rate(k - 1) (which tests one
      ! level too high and left lapse_rate_trop_level pointing at k_UM + 1).
      if (lapse_rate(k + 1) < lapse_trop .and. &
          lapse_rate(k)     > 0.0_r_def) then
        ! Lapse rate has dropped below the threshold. If this is maintained
        ! for 2km above then the WMO criteria for the tropopause has been
        ! met. Ported verbatim from UM's k2km loop in pws_tropoht_mod.F90:
        ! if the 2km-separated level would lie above heightcut_top, this
        ! candidate k is rejected outright (matching UM leaving tlev
        ! unset) rather than computing lapse_rate_above from a level
        ! above the search band. Likewise, if the column runs out before
        ! either bound is reached, the do loop simply ends without ever
        ! computing lapse_rate_above - UM has no "use whatever levels are
        ! left" fallback, so neither does this port.
        do kk = k + 1, nlayers
          if (height_wth(map_wth(1) + kk) > heightcut_top) exit
          dz = height_wth(map_wth(1) + kk) - height_wth(map_wth(1) + k)
          if (dz >= dz_trop) then
            lapse_rate_above = (t_wth(k) - t_wth(kk)) / dz
            if (lapse_rate_above < lapse_trop) then
              lapse_rate_trop_level = k
            end if
            exit
          end if
        end do
        if (lapse_rate_trop_level > 0) exit
      end if
    end if
  end do

  if (lapse_rate_trop_level > 0) then
    k = lapse_rate_trop_level

    ! Lapse rate for the interval below (k-1 ~ k) and above (k+1 ~ k+2) the
    ! tropopause level, used to interpolate height/temperature/pressure at
    ! the crossing point between the two lapse-rate lines.
    lapselwr = lapse_rate(k)
    lapseupr = (t_wth(k + 1) - t_wth(k + 2)) &
             / (height_wth(map_wth(1) + k + 2) - height_wth(map_wth(1) + k + 1))

    delta_lapse = lapselwr - lapseupr
    if (abs(delta_lapse) < vsmall) then
      if (delta_lapse >= 0.0_r_def) delta_lapse = vsmall
      if (delta_lapse <  0.0_r_def) delta_lapse = -vsmall
    end if

    trop_ht(map_2d(1)) =                                                   &
      ((t_wth(k)     + (lapselwr * height_wth(map_wth(1) + k)))   -        &
       (t_wth(k + 1) + (lapseupr * height_wth(map_wth(1) + k + 1)))) &
      / delta_lapse

    if (trop_ht(map_2d(1)) < height_wth(map_wth(1) + k)) then
      ! ensure trop height doesn't undershoot
      trop_ht(map_2d(1)) = height_wth(map_wth(1) + k)
    end if
    if (trop_ht(map_2d(1)) > height_wth(map_wth(1) + k + 1)) then
      ! or overshoot
      trop_ht(map_2d(1)) = height_wth(map_wth(1) + k + 1)
    end if

    if (.not. associated(trop_temp, empty_real_data)) then
      trop_temp(map_2d(1)) = t_wth(k) &
        - lapselwr * (trop_ht(map_2d(1)) - height_wth(map_wth(1) + k))
    end if

    ! Clamp lapselwr away from zero before it is used as a divisor below,
    ! unconditionally rather than nested inside the trop_press activation
    ! guard - matching UM's clamp (pws_tropoht_mod.F90, the ABS(lapselwr)
    ! < vsmall check following the tropopause-temperature calculation),
    ! which is likewise unconditional.
    if (abs(lapselwr) < vsmall) then
      if (lapselwr >= 0.0_r_def) lapselwr = vsmall
      if (lapselwr <  0.0_r_def) lapselwr = -vsmall
    end if

    if (.not. associated(trop_press, empty_real_data)) then
      ! Pressure at the tropopause is derived from the hydrostatic equation.
      press_at_k = p_zero * exner_in_wth(map_wth(1) + k)**(1.0_r_def / kappa)
      trop_press(map_2d(1)) = press_at_k &
        * (trop_temp(map_2d(1)) / t_wth(k))**(g_over_r / lapselwr)
    end if
  else
    ! No level in this column satisfies the WMO tropopause criteria (height
    ! + temperature band and sustained lapse-rate crossing). Matches UM's
    ! pws_tropoht, which leaves the diagnostic as missing data (tlev stays
    ! imdi) for that column rather than substituting a fallback level.
    trop_ht(map_2d(1)) = rmdi
    if (.not. associated(trop_temp, empty_real_data)) then
      trop_temp(map_2d(1)) = rmdi
    end if
    if (.not. associated(trop_press, empty_real_data)) then
      trop_press(map_2d(1)) = rmdi
    end if
  end if

  ! Reuses the shared ICAO standard-atmosphere conversion rather than
  ! duplicating it - see the @details note on icao_heights_kernel_code for
  ! why a kernel calling another module's procedure is acceptable here.
  ! icao_heights_kernel_code itself propagates trop_press == rmdi through
  ! to trop_icao_ht, so the missing-data case needs no separate handling.
  if (.not. associated(trop_icao_ht, empty_real_data)) then
    call icao_heights_kernel_code(nlayers, trop_icao_ht, trop_press, &
                                  g_over_r, ndf_2d, undf_2d, map_2d)
  end if

end subroutine tropopause_diags_code

end module tropopause_diags_kernel_mod
