! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module dry_adiabat_mod

implicit none

contains

! Subroutine to update the temperature of a parcel due to
! dry ascent or descent, under an adiabatic pressure change.
! Outputs the exner ratio to be used to rescale the temperature.


! Version for compressed arrays with all the condensed water
! species in a super-array
subroutine dry_adiabat( n_points, n_points_super,                              &
                        pressure_1, pressure_2,                                &
                        q_vap, q_cond_super, exner_ratio )

use comorph_constants_mod, only: real_cvprec, one, n_cond_species,             &
                                 l_approx_dry_adiabat, R_dry, R_vap, cp_dry
use set_cp_tot_mod, only: set_cp_tot

implicit none

! Number of points where we actually want to do something
integer, intent(in) :: n_points
! Number of points in the compressed super-array
! (to avoid having to repeatedly deallocate and reallocate when
!  re-using the same array, it maybe bigger than needed here).
integer, intent(in) :: n_points_super

! Pressure before
real(kind=real_cvprec), intent(in) :: pressure_1(n_points)
! Pressure after
real(kind=real_cvprec), intent(in) :: pressure_2(n_points)

! Water vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Super-array containing all condensed water mixing ratios
real(kind=real_cvprec), intent(in) :: q_cond_super                             &
                              ( n_points_super, n_cond_species )

! Exner ratio for rescaling temperature
real(kind=real_cvprec), intent(out) :: exner_ratio(n_points)

! Total heat capacity per unit dry-mass
real(kind=real_cvprec) :: cp_tot(n_points)

! Fractional change in pressure
real(kind=real_cvprec) :: dp_over_p(n_points)

! Gas constant over heat capacity per unit dry-mass
real(kind=real_cvprec) :: R_over_cp(n_points)

! Terms in Taylor expansion for (p2/p1)^(R/cp
real(kind=real_cvprec) :: term(n_points)
real(kind=real_cvprec) :: fac
real(kind=real_cvprec) :: n_real

! Loop counter
integer :: ic, n

! Number of terms to sum in Taylor expansion for (p2/p1)^(R/cp)
integer, parameter :: n_iter = 5
! Max fractional pressure change for the Taylor expansion to be accurate
real(kind=real_cvprec), parameter :: max_dp_over_p = 0.15_real_cvprec


! We expect that usually the fractional change in pressure will
! be small, in which case ( pressure_2 / pressure_1 )^R_over_cp
! can be accurately calculated using a Taylor expansion, which
! should be much cheaper than a call to power (**).
! So, for optimisation reasons, below we expand:
!
! (1 + dp/p)^(R/cp) = 1 + (R/cp) (dp/p) + (R/cp)(R/cp-1) 1/2 (dp/p)^2 + ...

do ic = 1, n_points
  ! Compute fractional change in pressure
  dp_over_p(ic) = pressure_2(ic) / pressure_1(ic) - one
  ! Initialise terms used in Taylor expansion...
  term(ic) = one
  exner_ratio(ic) = one
end do

if ( l_approx_dry_adiabat ) then
  ! Use approximate dry-adiabat, neglecting the effect of the
  ! water species on R and cp

  R_over_cp(1) = R_dry / cp_dry

  ! Taylor expansion:
  do n = 1, n_iter
    n_real = real(n,real_cvprec)
    fac = ( R_over_cp(1) - (n_real-one) ) / n_real
    do ic = 1, n_points
      term(ic) = term(ic) * fac * dp_over_p(ic)
      exner_ratio(ic) = exner_ratio(ic) + term(ic)
    end do
  end do

  ! Overwrite with call to power in the event that fractional pressure change
  ! is too large for the Taylor expansion to be accurate
  do ic = 1, n_points
    if ( abs(dp_over_p(ic)) > max_dp_over_p ) then
      exner_ratio(ic) = ( one + dp_over_p(ic) ) ** R_over_cp(1)
    end if
  end do

else
  ! Use full formula for the dry adiabat

  ! Calculate total parcel heat capacity per unit dry-mass
  call set_cp_tot( n_points, n_points_super, q_vap, q_cond_super,              &
                   cp_tot )

  do ic = 1, n_points
    R_over_cp(ic) = ( R_dry + q_vap(ic)*R_vap ) / cp_tot(ic)
  end do

  ! Taylor expansion:
  do n = 1, n_iter
    n_real = real(n,real_cvprec)
    do ic = 1, n_points
      fac = ( R_over_cp(ic) - (n_real-one) ) / n_real
      term(ic) = term(ic) * fac * dp_over_p(ic)
      exner_ratio(ic) = exner_ratio(ic) + term(ic)
    end do
  end do

  ! Overwrite with call to power in the event that fractional pressure change
  ! is too large for the Taylor expansion to be accurate
  do ic = 1, n_points
    if ( abs(dp_over_p(ic)) > max_dp_over_p ) then
      exner_ratio(ic) = ( one + dp_over_p(ic) ) ** R_over_cp(ic)
    end if
  end do

end if  ! ( l_approx_dry_adiabat )

return
end subroutine dry_adiabat



! Version for full 2-D arrays
subroutine dry_adiabat_2d( lb_p1, ub_p1, pressure_1, lb_p2, ub_p2, pressure_2, &
                           lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,               &
                           lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,               &
                           lb_s, ub_s, q_snow, lb_g, ub_g,q_graup,             &
                           exner_ratio )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full,                &
                     l_approx_dry_adiabat,                                     &
                     R_dry, R_vap, cp_dry
use set_cp_tot_mod, only: set_cp_tot_2d

implicit none

! Note: the array arguments below may have halos, and in here
! we don't care about the halos.
! We need to pass in the bounds of each array to use in
! the declarations, to ensure the do loops pick out the same
! indices of the arrays regardless of whether or not they
! have halos.  The lower and upper bounds are passed in through the
! argument list in the lb_* and ub_* integer arrays.  Those storing the
! bounds for 2D arrays must have 2 elements; one for each dimension
! of the array.

! Pressure before
integer, intent(in) :: lb_p1(2), ub_p1(2)
real(kind=real_hmprec), intent(in) :: pressure_1                               &
                                      ( lb_p1(1):ub_p1(1), lb_p1(2):ub_p1(2) )

! Pressure after
integer, intent(in) :: lb_p2(2), ub_p2(2)
real(kind=real_hmprec), intent(in) :: pressure_2                               &
                                      ( lb_p2(1):ub_p2(1), lb_p2(2):ub_p2(2) )

! Water species mixing ratios
integer, intent(in) :: lb_v(2), ub_v(2)
real(kind=real_hmprec), intent(in) :: q_vap                                    &
                                      ( lb_v(1):ub_v(1), lb_v(2):ub_v(2) )
integer, intent(in) :: lb_l(2), ub_l(2)
real(kind=real_hmprec), intent(in) :: q_cl                                     &
                                      ( lb_l(1):ub_l(1), lb_l(2):ub_l(2) )
integer, intent(in) :: lb_r(2), ub_r(2)
real(kind=real_hmprec), intent(in) :: q_rain                                   &
                                      ( lb_r(1):ub_r(1), lb_r(2):ub_r(2) )
integer, intent(in) :: lb_f(2), ub_f(2)
real(kind=real_hmprec), intent(in) :: q_cf                                     &
                                      ( lb_f(1):ub_f(1), lb_f(2):ub_f(2) )
integer, intent(in) :: lb_s(2), ub_s(2)
real(kind=real_hmprec), intent(in) :: q_snow                                   &
                                      ( lb_s(1):ub_s(1), lb_s(2):ub_s(2) )
integer, intent(in) :: lb_g(2), ub_g(2)
real(kind=real_hmprec), intent(in) :: q_graup                                  &
                                      ( lb_g(1):ub_g(1), lb_g(2):ub_g(2) )

! Exner ratio for rescaling temperature
real(kind=real_hmprec), intent(out) :: exner_ratio                             &
                                       ( nx_full, ny_full )

! Total heat capacity per unit dry-mass
real(kind=real_hmprec) :: cp_tot( nx_full, ny_full )

! Fractional change in pressure
real(kind=real_hmprec) :: dp_over_p( nx_full, ny_full )

! Gas constant over heat capacity per unit dry-mass
real(kind=real_hmprec) :: R_over_cp( nx_full, ny_full )

! Terms in Taylor expansion for (p2/p1)^(R/cp
real(kind=real_hmprec) :: term( nx_full, ny_full )
real(kind=real_hmprec) :: fac
real(kind=real_hmprec) :: n_real

! Copies of physical constants, converted to same precision as
! the full 2-D fields
real(kind=real_hmprec) :: R_dry_p, R_vap_p, cp_dry_p

! Loop counters
integer :: i, j, n

! Number of terms to sum in Taylor expansion for (p2/p1)^(R/cp)
integer, parameter :: n_iter = 10
! Max fractional pressure change for the Taylor expansion to be accurate
real(kind=real_hmprec), parameter :: max_dp_over_p = 0.15_real_hmprec
! One in native precision
real(kind=real_hmprec), parameter :: one = 1.0_real_hmprec

! Convert constants to same precision as the full 2-D fields
R_dry_p = real(R_dry,real_hmprec)
R_vap_p = real(R_vap,real_hmprec)
cp_dry_p = real(cp_dry,real_hmprec)

! We expect that usually the fractional change in pressure will
! be small, in which case ( pressure_2 / pressure_1 )^R_over_cp
! can be accurately calculated using a Taylor expansion, which
! should be much cheaper than a call to power (**).
! So, for optimisation reasons, below we expand:
!
! (1 + dp/p)^(R/cp) = 1 + (R/cp) (dp/p) + (R/cp)(R/cp-1) 1/2 (dp/p)^2 + ...

do j = 1, ny_full
  do i = 1, nx_full
    ! Compute fractional change in pressure
    dp_over_p(i,j) = pressure_2(i,j) / pressure_1(i,j) - one
    ! Initialise terms used in Taylor expansion...
    term(i,j) = one
    exner_ratio(i,j) = one
  end do
end do

if ( l_approx_dry_adiabat ) then
  ! Use approximate dry-adiabat, neglecting the effect of the
  ! water species on R and cp

  R_over_cp(1,1) = R_dry_p / cp_dry_p

  ! Taylor expansion:
  do n = 1, n_iter
    n_real = real(n,real_hmprec)
    fac = ( R_over_cp(1,1) - (n_real-one) ) / n_real
    do j = 1, ny_full
      do i = 1, nx_full
        term(i,j) = term(i,j) * fac * dp_over_p(i,j)
        exner_ratio(i,j) = exner_ratio(i,j) + term(i,j)
      end do
    end do
  end do

  ! Overwrite with call to power in the event that fractional pressure change
  ! is too large for the Taylor expansion to be accurate
  do j = 1, ny_full
    do i = 1, nx_full
      if ( abs(dp_over_p(i,j)) > max_dp_over_p ) then
        exner_ratio(i,j) = ( one + dp_over_p(i,j) ) ** R_over_cp(1,1)
      end if
    end do
  end do

else
  ! Use full formula for the dry adiabat

  ! Calculate total parcel heat capacity per unit dry-mass
  call set_cp_tot_2d( lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,                    &
                      lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,                    &
                      lb_s, ub_s, q_snow, lb_g, ub_g, q_graup,                 &
                      cp_tot )

  do j = 1, ny_full
    do i = 1, nx_full
      R_over_cp(i,j) = ( R_dry_p + q_vap(i,j)*R_vap_p ) / cp_tot(i,j)
    end do
  end do

  ! Taylor expansion:
  do n = 1, n_iter
    n_real = real(n,real_hmprec)
    do j = 1, ny_full
      do i = 1, nx_full
        fac = ( R_over_cp(i,j) - (n_real-one) ) / n_real
        term(i,j) = term(i,j) * fac * dp_over_p(i,j)
        exner_ratio(i,j) = exner_ratio(i,j) + term(i,j)
      end do
    end do
  end do

  ! Overwrite with call to power in the event that fractional pressure change
  ! is too large for the Taylor expansion to be accurate
  do j = 1, ny_full
    do i = 1, nx_full
      if ( abs(dp_over_p(i,j)) > max_dp_over_p ) then
        exner_ratio(i,j) = ( one + dp_over_p(i,j) ) ** R_over_cp(i,j)
      end if
    end do
  end do

end if  ! ( l_approx_dry_adiabat )

return
end subroutine dry_adiabat_2d


end module dry_adiabat_mod
