! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module solve_wf_cond_mod

implicit none

!----------------------------------------------------------------
! Indices to arguments in the super-array used to pass
! additional required arguments into brent_dekker_solve
!----------------------------------------------------------------

! Scalar inputs
integer, parameter :: i_rho_cond = 1
integer, parameter :: i_area_coef = 2
integer, parameter :: i_r_min_cond = 3

! Array inputs
integer, parameter :: i_q_cond = 1
integer, parameter :: i_n_cond = 2
integer, parameter :: i_rho_wet = 3
integer, parameter :: i_dt_over_lz = 4


contains


!----------------------------------------------------------------
! Subroutine to be called from inside Brent-Dekker routine to
! compute the error for each guess value of fall-speed wf_cond
!----------------------------------------------------------------
! Brent-Dekker will iteratively call this routine to converge on
! the value of wf_cond for which the output resid = 0.
! Note that this routine's argument list must follow the
! interface specified by f_of_x_interface in brent_dekker_mod.
subroutine resid_of_wf_cond( nc, n_points,                                     &
                             n_real_sca, n_real_arr,                           &
                             args_real_sca, args_real_arr,                     &
                             wf_cond, resid )

use comorph_constants_mod, only: real_cvprec

! Subroutines called inside here:
use fall_out_mod, only: fall_out
use set_cond_radius_mod, only: set_cond_radius
use fall_speed_mod, only: fall_speed

implicit none

! Number of points to do calculations at
integer, intent(in) :: nc

! Size of arrays
! (maybe larger than the actual number of points, as the
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

! Current guess value of the fall-speed
real(kind=real_cvprec), intent(in) :: wf_cond(n_points)

! Residual error (wf_after - wf_cond)
real(kind=real_cvprec), intent(out) :: resid(n_points)

! Loop counter
integer :: ic

! Work arrays
real(kind=real_cvprec) :: q_loc_cond(nc)
real(kind=real_cvprec) :: r_cond(nc)

! Mixing ratio after fall-out:
call fall_out( nc, wf_cond,                                                    &
               args_real_arr(:,i_dt_over_lz),                                  &
               args_real_arr(:,i_q_cond), q_loc_cond )
! Number conc. (constant so don't need to repeat)
! Particle radius:
call set_cond_radius( nc,                                                      &
                      args_real_sca(i_r_min_cond),                             &
                      args_real_sca(i_rho_cond),                               &
                      q_loc_cond,                                              &
                      args_real_arr(:,i_n_cond), r_cond )
! Fall-speed:
call fall_speed( nc,                                                           &
                 args_real_sca(i_area_coef),                                   &
                 args_real_sca(i_rho_cond),                                    &
                 args_real_arr(:,i_rho_wet),                                   &
                 r_cond, resid )
! (resulting fall-speed stored in resid)

! Calculate error wf_after - wf_cond
do ic = 1, n_points
  resid(ic) = resid(ic) - wf_cond(ic)
end do

return
end subroutine resid_of_wf_cond


!----------------------------------------------------------------
! Subroutine Sets up initial 2 guesses for the Brent-Decker
! algorithm, which will converge on a consistent fall-speed
!----------------------------------------------------------------
subroutine solve_wf_cond( n_points, nc, index_ic,                              &
                          cond_params, q_cond, n_cond,                         &
                          rho_wet, dt_over_lz,                                 &
                          q_loc_cond, r_cond, wf_cond, resid )

use comorph_constants_mod, only: real_cvprec, cond_params_type,                &
                                 solve_wf_n_iter, solve_wf_tolerance,          &
                                 i_solve_wf_check_converge, zero
use brent_dekker_mod, only: brent_dekker_solve
use fall_out_mod, only: fall_out
use set_cond_radius_mod, only: set_cond_radius

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where further iteration is actually needed
integer, intent(in) :: nc
! Indices of those points
integer, intent(in) :: index_ic(nc)

! Structure containing various fixed properties of the current
! condensed water species
type(cond_params_type), intent(in) :: cond_params

! Mixing ratio of hydrometeor before any fall-out
real(kind=real_cvprec), intent(in) :: q_cond(n_points)
! Number concentration before fall-out
real(kind=real_cvprec), intent(in) :: n_cond(n_points)
! Wet density of the air
real(kind=real_cvprec), intent(in) :: rho_wet(n_points)
! Factor delta_t / vert_len, used in precip fall-out calculation
real(kind=real_cvprec), intent(in) :: dt_over_lz(n_points)

! Mixing ratio after implicitly accounting for fall-out
real(kind=real_cvprec), intent(in out) :: q_loc_cond(n_points)
! Condensed water particle radius
real(kind=real_cvprec), intent(in out) :: r_cond(n_points)
! Fall-speed
real(kind=real_cvprec), intent(in out) :: wf_cond(n_points)
! Residual error of the fall-speed wf_after - wf_cond on input
real(kind=real_cvprec), intent(in) :: resid(n_points)
! Note: the final error is not output in this array


! Work arrays to pass into Brent-Dekker algorithm:

! First 2 guesses for wf_cond, and their corresponding errors
real(kind=real_cvprec) :: wf_cond_a(nc)
real(kind=real_cvprec) :: resid_a(nc)
real(kind=real_cvprec) :: wf_cond_b(nc)
real(kind=real_cvprec) :: resid_b(nc)

! Number of each sort of argument to pass in through Brent-Dekker
integer, parameter :: n_real_sca = 3  ! Real scalars (constants)
integer, parameter :: n_real_arr = 4  ! Real arrays

! Super-arrays to package the required array arguments together
real(kind=real_cvprec) :: args_real_sca(n_real_sca)
real(kind=real_cvprec) :: args_real_arr(nc,n_real_arr)

! Number of fields needed compressed onto work points
integer, parameter :: n_fields_cmpr = 7
! Addresses of those fields that aren't included in args_real_arr
integer, parameter :: i_q_loc_cond = 5
integer, parameter :: i_r_cond = 6
integer, parameter :: i_wf_cond = 7

! Super array to store all 7 fields compressed onto work points
! (only needed when iterations not done at all points)
real(kind=real_cvprec), allocatable :: super_cmpr(:,:)

! Loop counters
integer :: ic, ic2, i_field


! Package up the scalar arguments (constants) required
! by resid_of_wf_cond, called from inside brent_dekker_solve
args_real_sca(i_rho_cond)   = cond_params % rho
args_real_sca(i_area_coef)  = cond_params % area_coef
args_real_sca(i_r_min_cond) = cond_params % r_min

if ( nc == n_points ) then
  ! If iterations required at all points, work on the full arrays...

  ! Copy the initial guess passed in into guess b
  do ic = 1, n_points
    wf_cond_b(ic) = wf_cond(ic)
    resid_b(ic)   = resid(ic)
  end do

  ! Now, we know that the initial guess passed in will be an
  ! overestimate
  ! (it was calculated using the mixing-ratio without any
  !  fall-out accounted for; after fall-out, the mixing
  !  ratio will be smaller, leading to a lower fall-speed).
  ! To bracket the root, we need guess a to be an underestimate.
  ! An obvious, simple choice to go for is zero.  For a guess
  ! fall-speed of zero, there will be no fall-out, so that
  ! wf_after is equal to our first guess wf_cond.  Therefore
  ! the residual for this guess is wf_cond - zero:
  do ic = 1, n_points
    wf_cond_a(ic) = zero
    resid_a(ic)   = wf_cond(ic)
  end do

  ! Copy array arguments needed by resid_of_wf_cond, called
  ! from inside brent_dekker_solve
  do ic = 1, n_points
    args_real_arr(ic,i_q_cond)     = q_cond(ic)
    args_real_arr(ic,i_n_cond)     = n_cond(ic)
    args_real_arr(ic,i_rho_wet)    = rho_wet(ic)
    args_real_arr(ic,i_dt_over_lz) = dt_over_lz(ic)
  end do

  ! Call Brent-Dekker routine to converge on implicit solution
  ! for wf_cond
  call brent_dekker_solve( n_points, n_real_sca, n_real_arr,                   &
                           args_real_sca, args_real_arr,                       &
                           resid_of_wf_cond,                                   &
                           solve_wf_n_iter, solve_wf_tolerance,                &
                           i_solve_wf_check_converge,                          &
                           wf_cond_a, resid_a,                                 &
                           wf_cond_b, resid_b, wf_cond )

  ! Recalculate mixing ratio after fall-out and radius
  ! using final solved value of wf_cond
  call fall_out( n_points, wf_cond, dt_over_lz,                                &
                 q_cond, q_loc_cond )
  call set_cond_radius( n_points,                                              &
                        cond_params % r_min, cond_params % rho,                &
                        q_loc_cond, n_cond, r_cond )

else  ! ( nc == n_points )
  ! Not all points need iteration; do same as above but
  ! compress onto only points where iterations are needed

  allocate( super_cmpr( nc, n_fields_cmpr ) )

  do ic2 = 1, nc
    ic = index_ic(ic2)
    ! Compress inputs
    super_cmpr(ic2,i_q_cond)     = q_cond(ic)
    super_cmpr(ic2,i_n_cond)     = n_cond(ic)
    super_cmpr(ic2,i_rho_wet)    = rho_wet(ic)
    super_cmpr(ic2,i_dt_over_lz) = dt_over_lz(ic)
    ! Set first guesses for brent_dekker_solve
    wf_cond_b(ic2) = wf_cond(ic)
    resid_b(ic2)   = resid(ic)
    wf_cond_a(ic2) = zero
    resid_a(ic2)   = wf_cond(ic)
  end do

  ! Copy arguments input to brent_dekker_solve into temporary work array
  do i_field = 1, n_real_arr
    do ic2 = 1, nc
      args_real_arr(ic2,i_field) = super_cmpr(ic2,i_field)
    end do
  end do

  ! In call to Brent-Dekker, output new compressed solved value
  ! of wf_cond into work array super_cmpr
  call brent_dekker_solve( nc, n_real_sca, n_real_arr,                         &
                           args_real_sca, args_real_arr,                       &
                           resid_of_wf_cond,                                   &
                           solve_wf_n_iter, solve_wf_tolerance,                &
                           i_solve_wf_check_converge,                          &
                           wf_cond_a, resid_a,                                 &
                           wf_cond_b, resid_b, super_cmpr(:,i_wf_cond) )

  ! Recalculate mixing ratio after fall-out and radius
  ! using compressed work arrays resid_a, resid_b
  call fall_out( nc, super_cmpr(:,i_wf_cond), super_cmpr(:,i_dt_over_lz),      &
                 super_cmpr(:,i_q_cond), super_cmpr(:,i_q_loc_cond) )
  call set_cond_radius( nc,                                                    &
                        cond_params % r_min, cond_params % rho,                &
                        super_cmpr(:,i_q_loc_cond), super_cmpr(:,i_n_cond),    &
                        super_cmpr(:,i_r_cond) )

  ! Decompress solved outputs back to full arrays
  do ic2 = 1, nc
    ic = index_ic(ic2)
    wf_cond(ic)    = super_cmpr(ic2,i_wf_cond)
    q_loc_cond(ic) = super_cmpr(ic2,i_q_loc_cond)
    r_cond(ic)     = super_cmpr(ic2,i_r_cond)
  end do

  deallocate( super_cmpr )

end if  ! ( nc == n_points )


return
end subroutine solve_wf_cond


end module solve_wf_cond_mod
