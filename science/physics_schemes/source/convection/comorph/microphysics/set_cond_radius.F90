! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_cond_radius_mod

implicit none

contains


! Subroutine calculates the bulk radius of hydrometeors, given
! the mixing ratio and number concentration per unit dry-mass
! (assumes hydrometeors are spherical).
subroutine set_cond_radius( n_points, r_min_cond, rho_cond,                    &
                            q_cond, n_cond, r_cond )

use comorph_constants_mod, only: real_cvprec, sqrt_min_float,                  &
                                 pi, zero, one, third, four_thirds

implicit none

! Number of points
integer, intent(in) :: n_points

! Minimum particle radius for this condensed water species
real(kind=real_cvprec), intent(in) :: r_min_cond
! Density of the condensed water species
real(kind=real_cvprec), intent(in) :: rho_cond

! Mixing ratio of condensed water
real(kind=real_cvprec), intent(in) :: q_cond(n_points)
! Number concentration per unit dry-mass
real(kind=real_cvprec), intent(in) :: n_cond(n_points)

! Representative radius of hydrometeors
! (assuming they're all the same size and spherical!)
real(kind=real_cvprec), intent(out) :: r_cond(n_points)

! Scaling factors used to avoid floating-point underflow
real(kind=real_cvprec), parameter :: sqrt_min_float_recip                      &
                                     = one / sqrt_min_float
real(kind=real_cvprec), parameter :: sqrt_min_float_p1o3                       &
                                     = sqrt_min_float**third

! Loop counter
integer :: ic

do ic = 1, n_points
  ! Formula for radius of each particle:
  ! q = n rho 4/3 pi r^3
  ! => r = ( q / (n rho 4/3 pi) )^(1/3)
  !
  ! Note: the term inside the brackets (r^3) can get extremely small.
  ! In the most severe cases, it can lose precision or be truncated
  ! to zero due to floating-point underflow.
  ! Avoid this by redundantly scaling q_cond up by a factor and then
  ! dividing the result by factor^(1/3).
  r_cond(ic) = sqrt_min_float_p1o3 * ( (sqrt_min_float_recip * q_cond(ic))     &
               / ( four_thirds * pi * rho_cond * n_cond(ic) )                  &
               )**third
end do

! Apply a minimum value to r, consistent with size of CCN on
! which the hydrometeors grow
if ( r_min_cond > zero ) then
  do ic = 1, n_points
    r_cond(ic) = max( r_cond(ic), r_min_cond )
  end do
end if

return
end subroutine set_cond_radius

end module set_cond_radius_mod
