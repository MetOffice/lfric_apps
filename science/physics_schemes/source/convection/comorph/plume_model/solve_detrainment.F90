! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module solve_detrainment_mod

implicit none

! This module contains routines which are needed to solve the implicit
! detrainment formula, using an iterative root-finder
! (the Brent-Dekker algorithm)...

! We have an assumed linear function of buoyancy with x
! (where x is a dimensionless variable going from 0 to 1,
!  which just indicates where in the PDF we are;
!  0 is in the parcel core, 1 is in the existing parcel edge):
!
! Tv' = Tv'_core - (P+2)/(P+1) (Tv'_core - Tv'_mean) x
!
! The value of x where Tv' = 0 defines the new edge of the parcel x_edge;
! all air in the region x > x_edge is negatively buoyant and so
! is detrained.
! Rearranging the above with Tv'=0, we have:
!
! x_edge = (P+1)/(P+2) Tv'_core / (Tv'_core - Tv'_mean)
!
! Now the increase in environment Tv over the timestep due to
! compensating subsidence scales with the non-detrained fraction frac:
!
! dTv_env = delta_tv_k frac
!
! (the term delta_tv_k was precalculated in conv_level-step and stores
!  (Tv(k+1) - Tv(k) ) * mass-flux before detrainment / layer-mass).
!
! The buoyancies at end-of-timestep will be reduced by the compensating
! subsidence heating up the environment.  We seek an implicit-in-time
! solution such that the detrainment is consistent with the expected
! end-of-timestep buoyancies:
! Tv'_core_np1 = Tv'_core - dTv_env
! Tv'_mean_np1 = Tv'_mean - dTv_env
!
! Substituting these into the above expression for x_edge:
!
! x_edge = (P+1)/(P+2) (Tv'_core - delta_tv_k frac)
!                    / (Tv'_core - Tv'_mean)
!
! Now parcel mass-flux is assumed to be power-law distributed along x,
! so that the non-detrained fraction is given by:
!
! frac = x_edge^(P+1)
!
! Substituting this into the above formula for x_edge, we have:
!
! x_edge = (P+1)/(P+2) (Tv'_core - delta_tv_k x_edge^(P+1))
!                    / (Tv'_core - Tv'_mean)
!
! Annoyingly, there is no exact analytical rearrangement of this formula
! to calculate x_edge for a general PDF power P.
! (hence we're using an iterative root-finder algorithm to solve this).
!
! Rearranging the above, we can write:
!
! resid = (P+1)/(P+2) (Tv'_core - delta_tv_k x_edge^(P+1))
!                   / (Tv'_core - Tv'_mean)
!       - x_edge
!
! The correct solution will be the value of x_edge for which resid = 0.
!
! To solve this, we define a subroutine which calculates resid(x_edge),
! and pass that subroutine into the brent_dekker root-finder alorithm.


!----------------------------------------------------------------
! Indices to arguments in the super-array used to pass
! additional required arguments into brent_dekker_solve
!----------------------------------------------------------------

integer, parameter :: i_mean_buoy  = 1
integer, parameter :: i_core_buoy  = 2
integer, parameter :: i_power      = 3
integer, parameter :: i_delta_tv_k = 4


contains


!----------------------------------------------------------------
! Subroutine to be called from inside Brent-Dekker routine to
! compute the error for each guess value of x_edge
!----------------------------------------------------------------
! Brent-Dekker will iteratively call this routine to converge on
! the value of x_edge for which the output resid = 0.
! Note that this routine's argument list must follow the
! interface specified by f_of_x_interface in brent_dekker_mod.
subroutine resid_of_x_edge( nc, n_points,                                      &
                            n_real_sca, n_real_arr,                            &
                            args_real_sca, args_real_arr,                      &
                            x_edge, resid )

use comorph_constants_mod, only: real_cvprec, one

implicit none

! Number of points to do calculations at
integer, intent(in) :: nc

! Size of arrays
! (may be larger than the actual number of points, as the
!  same arrays are re-used even when some points have converged
!  and so have been removed from the list of points).
integer, intent(in) :: n_points

! Number of each sort of argument
integer, intent(in) :: n_real_sca   ! in real scalar arguments
integer, intent(in) :: n_real_arr   ! in real array arguments

! List of scalar inputs (constants)
real(kind=real_cvprec), intent(in) :: args_real_sca                            &
                                                ( n_real_sca )
! Super-array containing array arguments
real(kind=real_cvprec), intent(in) :: args_real_arr                            &
                                      ( n_points, n_real_arr )

! Current guess value of x_edge
real(kind=real_cvprec), intent(in) :: x_edge(n_points)

! Residual error
real(kind=real_cvprec), intent(out) :: resid(n_points)

! Store for power+1
real(kind=real_cvprec) :: power_p1

! Loop counter
integer :: ic

do ic = 1, nc
  power_p1 = args_real_arr(ic,i_power) + one
  ! Formula for the error resid as a function of x_edge
  ! (see the derivation at the top of this module):
  resid(ic) = (power_p1/(power_p1+one))                                        &
            * ( args_real_arr(ic,i_core_buoy)                                  &
                  - (x_edge(ic)**power_p1)                                     &
                    * args_real_arr(ic,i_delta_tv_k) )                         &
            / ( args_real_arr(ic,i_core_buoy)                                  &
              - args_real_arr(ic,i_mean_buoy) )                                &
            - x_edge(ic)
end do

return
end subroutine resid_of_x_edge


!----------------------------------------------------------------
! Subroutine Sets up initial 2 guesses for the Brent-Decker
! algorithm, which will converge on a consistent x_edge
!----------------------------------------------------------------
subroutine solve_detrainment( n_points, nc, index_ic,                          &
                              mean_buoy, core_buoy,                            &
                              power, delta_tv_k,                               &
                              x_edge, frac )

use comorph_constants_mod, only: real_cvprec, one,                             &
                                 solve_det_n_iter, solve_det_tolerance,        &
                                 i_solve_det_check_converge
use brent_dekker_mod, only: brent_dekker_solve

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where iterative solve is actually needed
integer, intent(in) :: nc
! Indices of those points
integer, intent(in) :: index_ic(nc)

! Mean buoyancy of the parcel
real(kind=real_cvprec), intent(in) :: mean_buoy(n_points)
! Core buoyancy of the parcel
real(kind=real_cvprec), intent(in) :: core_buoy(n_points)
! Power in the assumed power-law PDF
real(kind=real_cvprec), intent(in) :: power(n_points)
! Estimated change of environment virtual temperature at k
! due to compensating subsidence.
real(kind=real_cvprec), intent(in) :: delta_tv_k(n_points)

! Value of dimensionless variable x at the edge of the convection
! in: before the detrainment calculated here
! out: after the detrainment calculated here
real(kind=real_cvprec), intent(in out) :: x_edge(n_points)

! Non-detrained mass-fraction...
! in: before the detrainment calculated here
! out: after the detrainment calculated here
real(kind=real_cvprec), intent(in out) :: frac(n_points)

! First guess and its error
real(kind=real_cvprec) :: x_edge_a(nc)
real(kind=real_cvprec) :: resid_a(nc)

! Second guess and its error
real(kind=real_cvprec) :: x_edge_b(nc)
real(kind=real_cvprec) :: resid_b(nc)

! Compressed work copy of x_edge
real(kind=real_cvprec) :: x_edge_cmpr(nc)
! Compressed non-detrained fraction
real(kind=real_cvprec) :: frac_cmpr(nc)

! Number of each sort of argument to pass in through Brent-Dekker
integer, parameter :: n_real_sca = 0  ! Real scalars (constants)
integer, parameter :: n_real_arr = 4  ! Real arrays

! Super-arrays to package the required array arguments together
real(kind=real_cvprec) :: args_real_sca(n_real_sca)
real(kind=real_cvprec) :: args_real_arr(nc,n_real_arr)

! Note: there are no scalar arguments in this case, but the
! super-array for scalar args is still needed in the arg list.

! Loop counters
integer :: ic, ic2


! If iterations required at all points, just do a straight copy
if ( nc == n_points ) then

  ! Copy arguments into super-array
  do ic = 1, n_points
    args_real_arr(ic,i_mean_buoy)  = mean_buoy(ic)
    args_real_arr(ic,i_core_buoy)  = core_buoy(ic)
    args_real_arr(ic,i_power)      = power(ic)
    args_real_arr(ic,i_delta_tv_k) = delta_tv_k(ic)
  end do

  ! Set two initial guesses to bracket the root as closely as possible
  ! (subroutine declared at the bottom of this source file)
  call solve_det_init_guesses( n_points, x_edge, frac,                         &
                               mean_buoy, core_buoy, power, delta_tv_k,        &
                               x_edge_a, resid_a, x_edge_b, resid_b )

  ! Call  Brent-Dekker routine to converge on implicit solution
  ! for x_edge
  call brent_dekker_solve( n_points, n_real_sca, n_real_arr,                   &
                           args_real_sca, args_real_arr,                       &
                           resid_of_x_edge,                                    &
                           solve_det_n_iter, solve_det_tolerance,              &
                           i_solve_det_check_converge,                         &
                           x_edge_a, resid_a,                                  &
                           x_edge_b, resid_b, x_edge )

  ! Update the non-detrained fraction consistent with latest
  ! value of x_edge
  do ic = 1, n_points
    frac(ic) = x_edge(ic)**(power(ic)+one)
  end do

else  ! ( nc == n_points )
  ! Only detraining at some points...
  ! Do exactly as above but using the compression indices

  do ic2 = 1, nc
    ic = index_ic(ic2)
    args_real_arr(ic2,i_mean_buoy)  = mean_buoy(ic)
    args_real_arr(ic2,i_core_buoy)  = core_buoy(ic)
    args_real_arr(ic2,i_power)      = power(ic)
    args_real_arr(ic2,i_delta_tv_k) = delta_tv_k(ic)

    x_edge_cmpr(ic2) = x_edge(ic)
    frac_cmpr(ic2)   = frac(ic)
  end do

  call solve_det_init_guesses( nc, x_edge_cmpr, frac_cmpr,                     &
               args_real_arr(:,i_mean_buoy), args_real_arr(:,i_core_buoy),     &
               args_real_arr(:,i_power), args_real_arr(:,i_delta_tv_k),        &
               x_edge_a, resid_a, x_edge_b, resid_b )

  call brent_dekker_solve( nc, n_real_sca, n_real_arr,                         &
                           args_real_sca, args_real_arr,                       &
                           resid_of_x_edge,                                    &
                           solve_det_n_iter, solve_det_tolerance,              &
                           i_solve_det_check_converge,                         &
                           x_edge_a, resid_a,                                  &
                           x_edge_b, resid_b, x_edge_cmpr )

  do ic2 = 1, nc
    ic = index_ic(ic2)
    args_real_arr(ic2,i_power) = power(ic)
  end do
  do ic2 = 1, nc
    frac_cmpr(ic2) = x_edge_cmpr(ic2)**(args_real_arr(ic2,i_power)+one)
  end do
  do ic2 = 1, nc
    ic = index_ic(ic2)
    x_edge(ic) = x_edge_cmpr(ic2)
    frac(ic) = frac_cmpr(ic2)
  end do

end if  ! ( nc == n_points )


return
end subroutine solve_detrainment


! Routine to set the initial guesses for x_edge that bracket the root
! as closely as possible
subroutine solve_det_init_guesses( n_points, x_edge, frac,                     &
                                   mean_buoy, core_buoy, power, delta_tv_k,    &
                                   x_edge_a, resid_a, x_edge_b, resid_b )

use comorph_constants_mod, only: real_cvprec, zero, one, two,                  &
                                 sqrt_min_float, min_delta

implicit none

! Number of points
integer, intent(in) :: n_points

! Current-guess value of x_edge and the non-detrained fraction x_edge**(p+1)
real(kind=real_cvprec), intent(in) :: x_edge(n_points)
real(kind=real_cvprec), intent(in) :: frac(n_points)

! Mean buoyancy of the parcel
real(kind=real_cvprec), intent(in) :: mean_buoy(n_points)
! Core buoyancy of the parcel
real(kind=real_cvprec), intent(in) :: core_buoy(n_points)
! Power in the assumed power-law PDF
real(kind=real_cvprec), intent(in) :: power(n_points)
! Estimated change of environment virtual temperature at k
! due to compensating subsidence.
real(kind=real_cvprec), intent(in) :: delta_tv_k(n_points)

! First guess and its error
real(kind=real_cvprec), intent(out) :: x_edge_a(n_points)
real(kind=real_cvprec), intent(out) :: resid_a(n_points)

! Second guess and its error
real(kind=real_cvprec), intent(out) :: x_edge_b(n_points)
real(kind=real_cvprec), intent(out) :: resid_b(n_points)

! Store for precalculating p+1 / p+2
real(kind=real_cvprec) :: pp1_over_pp2(n_points)

! Temporary stores for alternative initial guess values of x_edge
real(kind=real_cvprec) :: x_edge_2, x_edge_3

! Loop counter
integer :: ic


! Precalculate common factor p+1 / p+2
do ic = 1, n_points
  pp1_over_pp2(ic) = (power(ic)+one) / (power(ic)+two)
end do

! Set 1st guess (a):
! Make this correspond to full detrainment (x_edge=0)
do ic = 1, n_points
  x_edge_a(ic) = zero
  ! Error is edge formula with frac set to zero, minus
  ! the guess which happens to be zero.
  resid_a(ic) = pp1_over_pp2(ic)                                               &
              * ( core_buoy(ic) )                                              &
              / ( core_buoy(ic) - mean_buoy(ic) )
end do

! Set 2nd guess (b):
! Make this correspond to some underestimate of the detrainment...
! Choose the highest detrainment rate (smallest x_edge) of 3 underestimates:
! 1) The current input value of x_edge
! 2) The explicit estimate of x_edge (ignoring the delta_tv_k term)
! 3) The limit as mean_buoy => core_buoy (only the delta_tv_k term)
do ic = 1, n_points
  ! Initialise using input value of x_edge and frac
  x_edge_b(ic) = x_edge(ic)
  ! Set error = edge formula minus value of edge used
  resid_b(ic) = pp1_over_pp2(ic)                                               &
              * ( core_buoy(ic) - frac(ic) * delta_tv_k(ic) )                  &
              / ( core_buoy(ic) - mean_buoy(ic) )                              &
              - x_edge(ic)
end do

do ic = 1, n_points
  ! Guess 2) explicit estimate
  !          (+ min_delta avoids failing to bracket root due to rounding error)
  x_edge_2 = pp1_over_pp2(ic)                                                  &
           * core_buoy(ic) / ( core_buoy(ic) - mean_buoy(ic) ) + min_delta
  ! Guess 3) limit as mean_buoy => core_buoy
  x_edge_3 = ( core_buoy(ic) / max( delta_tv_k(ic), sqrt_min_float )           &
             )**(one/(power(ic)+one))
  if ( x_edge_2 < x_edge(ic) .and. x_edge_2 <= x_edge_3 ) then
    ! Use guess 2
    x_edge_b(ic) = x_edge_2
    ! Set error = edge formula minus value of edge used
    resid_b(ic) = pp1_over_pp2(ic)                                             &
                * ( core_buoy(ic) - ( x_edge_2**(power(ic)+one) )              &
                                    * delta_tv_k(ic) )                         &
                / ( core_buoy(ic) - mean_buoy(ic) )                            &
                - x_edge_2
  else if ( x_edge_3 < x_edge(ic) ) then
    ! Use guess 3
    x_edge_b(ic) = x_edge_3
    ! Set error; substituting frac = x_edge**(p+1) = core_buoy/delta_tv_k,
    ! note that the numerator in the edge formula cancels out and goes to
    ! zero, so the error is just -x_edge.
    resid_b(ic) = -x_edge_3
  end if
end do


return
end subroutine solve_det_init_guesses


end module solve_detrainment_mod
