! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module fall_speed_mod

implicit none

contains


! Subroutine calculates the fall-speed of a hydrometeor with
! a specified radius, density and drag-coefficient.
subroutine fall_speed( n_points, area_coef, rho_cond, rho_air,                 &
                       r_cond, wf_cond )

use comorph_constants_mod, only: real_cvprec, sqrt_min_delta,                  &
                                 gravity, kin_visc, drag_coef_cond,            &
                                 zero, one, two, six, half, four_thirds
                                 ! (nope; this subroutine can't count to 3)

implicit none

! Number of points
integer, intent(in) :: n_points

! Ratio of hydrometeor actual cross-section area to what you get
! by assuming a sphere.  This scales the drag force, making it
! bigger for non-spherical species such as ice crystals
real(kind=real_cvprec), intent(in) :: area_coef

! Density of the condensed water particles
real(kind=real_cvprec), intent(in) :: rho_cond

! Density of the surround ing air
real(kind=real_cvprec), intent(in) :: rho_air(n_points)

! Radius of the condensed water particles
real(kind=real_cvprec), intent(in) :: r_cond(n_points)

! Output fall-speed
real(kind=real_cvprec), intent(out) :: wf_cond(n_points)

! Precalculate some common non-array factors
real(kind=real_cvprec), parameter :: six_kin_visc = six * kin_visc
real(kind=real_cvprec), parameter :: six_kin_visc_sq = six_kin_visc            &
                                                     * six_kin_visc

! Loop counter
integer :: ic


! We assume the drag acting on a hydrometeor is just the sum of
! Stoke's Law (the low Reynolds number limit) and the quadratic drag
! equation (the high Reynolds number limit), for a sphere.
! The sum of these two drag forces must balance acceleration due to buoyancy:
!
! Fd/m  =  acf (rho_a/rho_c) ( 9/2 nu wf / R^2  +  3/8 Cd wf^2 / R )
!       =  g (rho_c - rho_a)/rho_c
!
! rho_a = air density
! rho_c = hydrometeor density
! acf = dimensionless factor accounting for non-spherical shape
! nu = molecular viscosity
! Cd = dimensionless drag coefficient
! R = hydrometeor radius
! wf = fall-speed
!
! This rearranges into the following quadratic equation for the fall-speed wf:
!
! 1/2 wf^2  +  6 nu / (Cd R) wf  -  4/3 R / (Cd acf) g (rho_c/rho_a - 1)  =  0
!
! This has the solution:
! wf = -b + sqrt( b^2 - 2c )
!    = -6 nu / (Cd R)  +  sqrt( ( 6 nu / (Cd R) )^2
!                             + 2 4/3 R / (Cd acf) g (rho_c/rho_a - 1) )
!    = 6 nu / (Cd R) ( -1 + sqrt( 1 + 2 Cd g (rho_c/rho_a - 1)
!                                     1/( (6 nu)^2 acf ) 4/3 R^3 ) )

! For small Reynolds number, the full formula is unsafe
! (result is small residual of large opposing terms).

! We need to compute something of the form
!   y = sqrt( 1 + x ) - 1
! which is a small residual between compensating terms when x is
! small, and so loses precision.  When x is small, the Taylor
! expansion can be used: y = 1/2 x.
! We use a threshold in x for whether to use full formula
! or Taylor.  We try to set this near the value of x where the
! precision of y from the full formula falls below that from the
! Taylor expansion.  Used threshold is sqrt(epsilon).

do ic = 1, n_points
  ! Calculate the term x first and see if it is big enough
  ! to yield a safely calculable residual:
  wf_cond(ic) = two * drag_coef_cond * gravity                                 &
                * ( rho_cond/rho_air(ic) - one )                               &
                * ( one / ( six_kin_visc_sq * area_coef ) )                    &
                * four_thirds * r_cond(ic)*r_cond(ic)*r_cond(ic)
end do

do ic = 1, n_points
  if ( wf_cond(ic) > sqrt_min_delta ) then
     ! The term x is large enough; use full formula
    wf_cond(ic) = sqrt( one + wf_cond(ic) ) - one
  else
     ! The term x is too small; use Taylor expansion
    wf_cond(ic) = half * wf_cond(ic)
  end if
end do

do ic = 1, n_points
  ! Divide by particle radius (taking care in case r_cond = 0)
  if ( r_cond(ic) > zero ) wf_cond(ic) = wf_cond(ic) / r_cond(ic)
  ! Note that if r_cond is exactly zero, wf_cond will also
  ! already be zero as it has a factor of r_cond**3 in it above.
end do

do ic = 1, n_points
  ! Complete the formula to obtain fall-speed
  wf_cond(ic) = ( six_kin_visc / drag_coef_cond ) * wf_cond(ic)
end do


return
end subroutine fall_speed

end module fall_speed_mod
