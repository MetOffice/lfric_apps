! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module wind_w_eqn_mod

implicit none

contains

! Subroutine integrates the parcel vertical velocity equation
subroutine wind_w_eqn( n_points,                                               &
                       height_prev, height_1, height_2,                        &
                       tv_ex_1, tv_ex_2, delta_tv_1, delta_tv_2,               &
                       env_tv_1, env_tv_2,                                     &
                       par_w_drag, buoydz, wind_w_ex )

use comorph_constants_mod, only: real_cvprec, sqrt_min_float, half, two,       &
                                 gravity

implicit none

! Number of points
integer, intent(in) :: n_points

! Height of start of level-step, and latest two sub-levels
real(kind=real_cvprec), intent(in) :: height_prev(n_points)
real(kind=real_cvprec), intent(in) :: height_1(n_points)
real(kind=real_cvprec), intent(in) :: height_2(n_points)

! Parcel virtual temperature excess on latest two sub-levels
real(kind=real_cvprec), intent(in) :: tv_ex_1(n_points)
real(kind=real_cvprec), intent(in) :: tv_ex_2(n_points)

! Environment Tv increments due to subsidence term
real(kind=real_cvprec), intent(in) :: delta_tv_1(n_points)
real(kind=real_cvprec), intent(in) :: delta_tv_2(n_points)

! Environment Tv
real(kind=real_cvprec), intent(in) :: env_tv_1(n_points)
real(kind=real_cvprec), intent(in) :: env_tv_2(n_points)

! Drag coefficient / s-1
real(kind=real_cvprec), intent(in) :: par_w_drag(n_points)

! Buoyancy integrated over sub-levels (updated for current sub-level here)
real(kind=real_cvprec), intent(in out) :: buoydz(n_points)

! Parcel vertical velocity excess (updated with buoyancy and drag here)
real(kind=real_cvprec), intent(in out) :: wind_w_ex(n_points)

! Buoyancies in ms-2
real(kind=real_cvprec) :: buoy_1
real(kind=real_cvprec) :: buoy_2

! Increment to w' due to drag
real(kind=real_cvprec) :: drag_inc

! Work variables storing w' and 1/2 w'^2
real(kind=real_cvprec) :: wex
real(kind=real_cvprec) :: vert_ke

! Terminal velocity of the parcel under current drag and buoydz
real(kind=real_cvprec) :: wt

! Loop counter
integer :: ic


! We have:
! Dwp/Dt = b' - d ( wp - we )
!
! where b' is buoyancy / ms-2, and d is the drag coefficient / s-1.
!
! Now we take a frame of reference following the environment
! (assuming it is defined in a Lagrangian sense, following we):
!
! D/Dt = d/dt + (w-we) d/dz
!
! So, following the parcel,
! Dwe/Dt = (wp-we) dwe/dz
! Dwp/Dt = (wp-we) dwp/dz
! (assumed to be in an Eulerian steady state, so d/dt=0)
!
! Putting it all together:
!
! D/Dt(wp-we) = (wp-we) d/dz(wp-we)
!             = Dwp/Dt - (wp-we) dwe/dz
!             = b' - d (wp-we) - (wp-we) dwe/dz
!
! => d/dz(wp-we) = 1/(wp-we) b' - d - dwe/dz
!
! Now, the first term on the r.h.s. is easiest to integrate by rearranging
! the equation in terms of w'^2:
! d/dz( 1/2 w'^2 ) = b' - w' ( d + dwe/dz )
!
! But the 2nd term on the r.h.s. naturally just gives a linear change
! in w' with height.
!
! Therefore we discretize the equation by integrating buoyancy
! to increment w'^2, and then converting back to w' to subtract
! the drag and dwe/dz terms.
!
! Note that the dwe/dz term is accounted for by the fact that
! we recalculated w' at the next level-step using the environment w
! from the next level, so we don't need to explicitly add it on.


do ic = 1, n_points

  ! Compute buoyancy, accounting for change in env Tv due to subsidence
  buoy_1 = gravity * ( tv_ex_1(ic) - delta_tv_1(ic) )                          &
                   / max( env_tv_1(ic), sqrt_min_float )
  buoy_2 = gravity * ( tv_ex_2(ic) - delta_tv_2(ic) )                          &
                   / max( env_tv_2(ic), sqrt_min_float )

  ! Integrate buoyancy up to current sub-level
  buoydz(ic) = buoydz(ic)                                                      &
             + half * ( buoy_1 + buoy_2 ) * abs( height_2(ic) - height_1(ic) )
  ! Want this to have same sign as the buoyancies, including for
  ! downdrafts where dz < 0, so take ABS of dz

  ! Integrate drag
  drag_inc = -par_w_drag(ic) * ( height_2(ic) - height_prev(ic) )
  ! Note that the drag coefficient is always guaranteed to be positive.
  ! It is then scaled by -dz, so will always be a negative increment
  ! for updrafts, positive increment for downdrafts.

  ! Add on half the drag integrated to current sub-level height
  wex = wind_w_ex(ic) + half*drag_inc

  ! Convert to (signed) KE and add on buoyancy term
  ! (negative buoyancy should increase negative KE for downdrafts)
  vert_ke = half*wex*abs(wex) + buoydz(ic)

  ! Convert back to w-excess (retain the sign if KE goes negative),
  ! and add the other half of the drag
  wex = sign( sqrt(abs(two*vert_ke)), vert_ke ) + half*drag_inc

  ! w' should converge towards buoydz / (-drag_inc).
  ! Discretisation errors can sometimes cause it to over- or under-
  ! shoot this value and oscillate.  Avoid this by not allowing the
  ! buoyancy+drag increment to make it cross from under to over this value
  ! (or vice-versa)
  wt = buoydz(ic) / sign( max( abs(drag_inc), sqrt_min_float ), -drag_inc )
  if ( wind_w_ex(ic) > wt ) then
    wex = max( wex, wt )
  else
    wex = min( wex, wt )
  end if

  ! Output final value
  wind_w_ex(ic) = wex

end do


return
end subroutine wind_w_eqn

end module wind_w_eqn_mod
