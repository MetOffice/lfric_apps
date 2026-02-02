! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module init_sublevs_mod

implicit none

contains

! Subroutine to initialise the values of various fields held
! at multiple sub-level steps in a super-array in conv_level-step
subroutine init_sublevs( n_points, n_points_super, l_down,                     &
                         prev_height, next_height,                             &
                         env_prev_super, env_next_super,                       &
                         par_mean_virt_temp, par_core_virt_temp,               &
                         par_mean_wind_w, par_core_wind_w,                     &
                         par_prev_super, par_conv_super,                       &
                         delta_tv, sum_massflux_det,                           &
                         sublevs, i_next, i_sat, i_core_sat )

use comorph_constants_mod, only: real_cvprec, zero, min_float,                 &
                                 l_par_core
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_height, j_massflux_d, j_env_tv, j_delta_tv, j_env_w,  &
                       j_mean_buoy, j_core_buoy,                               &
                       j_mean_wex, j_core_wex

use parcel_type_mod, only: n_par, i_massflux_d
use env_half_mod, only: n_env_half, i_virt_temp, i_wind_w_half

implicit none

! Number of points
integer, intent(in) :: n_points

! Array size of some super-arrays
integer, intent(in) :: n_points_super

! Flag for downdrafts vs updrafts
logical, intent(in) :: l_down

! Height at start and end of the level-step
real(kind=real_cvprec), intent(in) :: prev_height(n_points)
real(kind=real_cvprec), intent(in) :: next_height(n_points)

! Super-arrays containing any environment fields required
! at start and end of the current half-level step
real(kind=real_cvprec), intent(in) :: env_prev_super                           &
                                      ( n_points_super, n_env_half )
real(kind=real_cvprec), intent(in) :: env_next_super                           &
                                      ( n_points_super, n_env_half )

! Parcel mean and core virtual temperature at start of the level-step
real(kind=real_cvprec), intent(in) :: par_mean_virt_temp(n_points)
real(kind=real_cvprec), intent(in) :: par_core_virt_temp(n_points)

! Parcel mean and core vertical velocity at start of the level-step
real(kind=real_cvprec), intent(in) :: par_mean_wind_w(n_points)
real(kind=real_cvprec), intent(in) :: par_core_wind_w(n_points)

! Parcel properties at start and end of the level-step
real(kind=real_cvprec), intent(in) :: par_prev_super ( n_points, n_par )
real(kind=real_cvprec), intent(in) :: par_conv_super ( n_points_super, n_par )

! Difference in virtual temperature between level k and next full-level,
! divided by layer-mass (used to pre-estimate subsidence increment to Tv)
real(kind=real_cvprec), intent(in) :: delta_tv( n_points )

! Effective mass-flux at start of this half-level step for detrainment
real(kind=real_cvprec), intent(in) :: sum_massflux_det(n_points)

! Super-array storing the parcel buoyancies and other properties
! at up to 4 sub-level heights within the current level-step:
! a) Start of the level-step
! b) End of the level-step
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(out) :: sublevs                                 &
                                       ( n_points, n_sublev_vars, max_sublevs )

! Address of next model-level interface in sublevs
integer, intent(out) :: i_next(n_points)
! Address of saturation height in sublevs
integer, intent(out) :: i_sat(n_points)
integer, intent(out) :: i_core_sat(n_points)

! Loop counters
integer :: ic, i_field, i_lev


do ic = 1, n_points

  ! Initialise sub-level height addresses
  ! (no intervening saturation height levels found yet)
  i_next(ic) = i_prev + 1
  i_sat(ic) = 0
  i_core_sat(ic) = 0

  ! Initialise heights, mass-fluxes, and env Tv and w for start and end of
  ! this level-step.  Note the "next" mass-flux so far only accounts for
  ! entrainment; detrainment will be accounted for in set_det.
  sublevs(ic,j_height,i_prev)         = prev_height(ic)
  sublevs(ic,j_height,i_next(ic))     = next_height(ic)
  sublevs(ic,j_massflux_d,i_prev)     = par_prev_super(ic,i_massflux_d)
  sublevs(ic,j_massflux_d,i_next(ic)) = par_conv_super(ic,i_massflux_d)
  sublevs(ic,j_env_tv,i_prev)         = env_prev_super(ic,i_virt_temp)
  sublevs(ic,j_env_tv,i_next(ic))     = env_next_super(ic,i_virt_temp)
  sublevs(ic,j_env_w,i_prev)          = env_prev_super(ic,i_wind_w_half)
  sublevs(ic,j_env_w,i_next(ic))      = env_next_super(ic,i_wind_w_half)

  ! Initialise parcel mean buoyancy and w-excess at start of the level-step
  ! (values at next will be set in parcel_dyn; initialise to zero for now)
  sublevs(ic,j_mean_buoy,i_prev) = par_mean_virt_temp(ic)                      &
                                 - env_prev_super(ic,i_virt_temp)
  sublevs(ic,j_mean_wex,i_prev)  = par_mean_wind_w(ic)                         &
                                 - env_prev_super(ic,i_wind_w_half)
  sublevs(ic,j_mean_buoy,i_next(ic)) = zero
  sublevs(ic,j_mean_wex,i_next(ic))  = zero

end do

if ( l_par_core ) then
  ! Initialise core buoyancy and w-excess at start of the level-step, if used
  do ic = 1, n_points
    sublevs(ic,j_core_buoy,i_prev) = par_core_virt_temp(ic)                    &
                                   - env_prev_super(ic,i_virt_temp)
    sublevs(ic,j_core_wex,i_prev)  = par_core_wind_w(ic)                       &
                                   - env_prev_super(ic,i_wind_w_half)
    sublevs(ic,j_core_buoy,i_next(ic)) = zero
    sublevs(ic,j_core_wex,i_next(ic))  = zero
  end do
end if

! Estimate change in environment virtual temperature
! due to compensating subsidence...
do ic = 1, n_points

  ! Now compute the new value at the current level-step.
  ! This is the subsided value of Tv from the next full-level
  ! minus the value at level k,
  ! * dt * massflux / layer-mass
  ! We also scale by the implicitness weight alpha_detrain.
  ! (the input delta_tv stores the Tv difference divided by layer-mass,
  !  scaled by dt * alpha_detrain)
  ! Note we want the mass-flux at next, whereas sum_massflux_det is at prev,
  ! so scaling up by ratio next/prev mass-fluxes to account for entrainment.
  sublevs(ic,j_delta_tv,i_next(ic)) = delta_tv(ic)                             &
                      * sum_massflux_det(ic)                                   &
                      * ( par_conv_super(ic,i_massflux_d)                      &
                   / max( par_prev_super(ic,i_massflux_d), min_float ) )

end do

! For safety, reset delta_tv_next to zero in statically-unstable
! layers so we just use explicit value of env Tv for detrainment.
! Expected sign depends on whether this is updraft or downdraft
if ( l_down ) then
  do ic = 1, n_points
    sublevs(ic,j_delta_tv,i_next(ic))                                          &
      = min( sublevs(ic,j_delta_tv,i_next(ic)), zero )
  end do
else
  do ic = 1, n_points
    sublevs(ic,j_delta_tv,i_next(ic))                                          &
      = max( sublevs(ic,j_delta_tv,i_next(ic)), zero )
  end do
end if

! TEMPORARY CODE TO PRESERVE KGO; TO BE REMOVED SOON:
! Set delta_tv_prev equal to value at next (assume delta_tv const with height)
do ic = 1, n_points
  sublevs(ic,j_delta_tv,i_prev) = sublevs(ic,j_delta_tv,i_next(ic))
end do

! The estimated environment Tv after compensating subsidence
! can now be calculated as
!   virt_temp_impl = env_virt_temp + delta_tv * frac
!
! where frac is the non-detrained fraction of the mass-flux.
!
! Note that this estimate of the final Tv does NOT
! account for the increments from:
! - entrainment and detrainment.
! - other convective drafts.
! - boundary-layer mixing...
! It is therefore not expected to yield an accurate implicit
! solution.  It is only designed to be a relatively simple
! modification of the explicit detrainment calculation, to
! introduce just enough implicitness to smooth out unphysical
! noisy / intermittent behaviour from one timestep to the next.

! Initialise sub-level data in not-yet-used sub-level height addresses to zero
! (most of the time not all the sub-levels are used, but some calculations
!  are best written to loop over all of them; setting the data at unused
!  sub-levels to zero ensures this is safe!)
do i_lev = i_prev+2, max_sublevs
  do i_field = 1, n_sublev_vars
    do ic = 1, n_points
      sublevs(ic,i_field,i_lev) = zero
    end do
  end do
end do


return
end subroutine init_sublevs

end module init_sublevs_mod
