! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_cond_properties_mod

implicit none

contains


! Subroutine solves the mixing-ratio / number concentration
! relationship to find number conc. n and particle radius r,
! calculates the fall-speed wf, and makes an initial estimate
! for the mixing ratio (and consistent n, r) after fall-out
! of the hydrometeor.
! An implicit solution is found such that n, r, and wf (and
! the amount of fall-out occuring with speed wf) are
! consistent with the hydrometeor's mixing ratio after fall-out.
!
! The particle radius, number concentration and fall-speed
! are then used to calculate coefficients for the exchange
! of water vapour and heat between the condensed water and
! the surrounding air
subroutine calc_cond_properties( n_points,                                     &
                                 cond_params, ref_temp, q_cond,                &
                                 rho_dry, rho_wet, dt_over_lz,                 &
                                 q_loc_cond, n_cond,                           &
                                 r_cond, wf_cond,                              &
                                 kq_cond, kt_cond )

use comorph_constants_mod, only: real_cvprec, cond_params_type,                &
                     melt_temp, homnuc_temp, fac_tdep_n,                       &
                     solve_wf_tolerance
use set_cond_radius_mod, only: set_cond_radius
use fall_speed_mod, only: fall_speed
use fall_out_mod, only: fall_out
use solve_wf_cond_mod, only: solve_wf_cond
use calc_kqkt_mod, only: calc_kqkt

implicit none

! Number of points
integer, intent(in) :: n_points

! Structure containing various fixed properties of the current
! condensed water species
type(cond_params_type), intent(in) :: cond_params

! Reference temperature
real(kind=real_cvprec), intent(in) :: ref_temp(n_points)

! Mixing ratio of hydrometeor before any fall-out
real(kind=real_cvprec), intent(in) :: q_cond(n_points)

! Dry-density and wet-density of the air
real(kind=real_cvprec), intent(in) :: rho_dry(n_points)
real(kind=real_cvprec), intent(in) :: rho_wet(n_points)

! Factor delta_t / vert_len, used in precip fall-out calculation
real(kind=real_cvprec), intent(in) :: dt_over_lz(n_points)

! Mixing ratio of hydrometeor after an initial guess for the
! amount of fall-out.  The particle radius and fall-speed are
! set consistent with this value
real(kind=real_cvprec), intent(out) :: q_loc_cond(n_points)
! Number concentration of particles
real(kind=real_cvprec), intent(out) :: n_cond(n_points)
! Condensed water particle radius
real(kind=real_cvprec), intent(out) :: r_cond(n_points)
! Fall-speed
real(kind=real_cvprec), intent(out) :: wf_cond(n_points)

! Coeficients for exchange of water vapour and heat between
! the condensed water particles and the surrounding air
real(kind=real_cvprec), intent(out) :: kq_cond(n_points)
real(kind=real_cvprec), intent(out) :: kt_cond(n_points)

! Error of first-guess fall-speed
real(kind=real_cvprec) :: resid(n_points)

! Number of points where iteration is required to converge on
! implicit fall-speed solution
integer :: nc
! Indices of those points
integer :: index_ic(n_points)

! Loop counter
integer :: ic


! Set the first-guess number concentration
do ic = 1, n_points
  ! Currently just set to a constant value for each species
  n_cond(ic) = cond_params % n
end do

! If this species uses temperature-dependent number concentration
if ( cond_params % l_tdep_n ) then
  do ic = 1, n_points
    ! Scale n by temperature-dependent factor
    n_cond(ic) = n_cond(ic) * exp( fac_tdep_n                                  &
       * ( max(min( ref_temp(ic), melt_temp ),homnuc_temp)                     &
         - melt_temp )  )
  end do
end if

! Set first-guess particle radius
call set_cond_radius( n_points,                                                &
                      cond_params % r_min, cond_params % rho,                  &
                      q_cond, n_cond, r_cond )

! Set first-guess fall-speed
call fall_speed( n_points,                                                     &
                cond_params % area_coef, cond_params % rho,                    &
                rho_wet, r_cond, wf_cond )


! Find the fall-speed that you'd get by using the current-guess
! fall-speed to predict fall-out, final mixng ratio, and
! particle radius consistent with that:

! Mixing ratio after fall-out:
call fall_out( n_points, wf_cond, dt_over_lz,                                  &
               q_cond, q_loc_cond )
! Number conc. (constant so don't need to repeat)
! Particle radius:
call set_cond_radius( n_points,                                                &
                      cond_params % r_min, cond_params % rho,                  &
                      q_loc_cond, n_cond, r_cond )
! Fall-speed:
call fall_speed( n_points,                                                     &
                 cond_params % area_coef, cond_params % rho,                   &
                 rho_wet, r_cond, resid )
! (resulting fall-speed stored in resid)

! Find the error in the implicit fall-speed = final wf - guess wf
do ic = 1, n_points
  resid(ic) = resid(ic) - wf_cond(ic)
end do


! Set points where further iterations are needed to converge on
! resid = wf_after - wf_cond = 0.
nc = 0
do ic = 1, n_points
  ! wf_cond needs to be implicitly corrected at points where the
  ! fractional error is greater than a set tolerance:
  if ( abs(resid(ic)) > solve_wf_tolerance*wf_cond(ic) ) then
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do

if ( nc > 0 ) then

  call solve_wf_cond( n_points, nc, index_ic,                                  &
                      cond_params, q_cond, n_cond,                             &
                      rho_wet, dt_over_lz,                                     &
                      q_loc_cond, r_cond, wf_cond, resid )

end if


! Calculate coefficients for exchange of water vapour and heat
! between the condensed water and the surrounding air
call calc_kqkt( n_points, cond_params % area_coef,                             &
                n_cond, r_cond, wf_cond, rho_dry,                              &
                kq_cond, kt_cond )


return
end subroutine calc_cond_properties



! Compression wrapper for the above
subroutine calc_cond_properties_cmpr( n_points, nc, index_ic,                  &
                                      cond_params,                             &
                                      ref_temp, q_cond,                        &
                                      rho_dry, rho_wet,                        &
                                      dt_over_lz,                              &
                                      q_loc_cond, n_cond,                      &
                                      r_cond, wf_cond,                         &
                                      kq_cond, kt_cond )

use comorph_constants_mod, only: real_cvprec, cond_params_type, zero,          &
                                 cmpr_thresh

implicit none

! Number of points in full arrays
integer, intent(in) :: n_points
! Number of points where condensed water field is non-zero
integer, intent(in) :: nc
! Indices of non-zero points
integer, intent(in) :: index_ic(nc)

! Other arguments identical to calc_cond_properties above
type(cond_params_type), intent(in) :: cond_params
real(kind=real_cvprec), intent(in) :: ref_temp(n_points)
real(kind=real_cvprec), intent(in) :: q_cond(n_points)
real(kind=real_cvprec), intent(in) :: rho_dry(n_points)
real(kind=real_cvprec), intent(in) :: rho_wet(n_points)
real(kind=real_cvprec), intent(in) :: dt_over_lz(n_points)
real(kind=real_cvprec), intent(in out) :: q_loc_cond(n_points)
real(kind=real_cvprec), intent(in out) :: n_cond(n_points)
real(kind=real_cvprec), intent(in out) :: r_cond(n_points)
real(kind=real_cvprec), intent(in out) :: wf_cond(n_points)
real(kind=real_cvprec), intent(in out) :: kq_cond(n_points)
real(kind=real_cvprec), intent(in out) :: kt_cond(n_points)
! Outputs need intent inout so that values initialised to zero
! are preserved at points that aren't in the compression list.

! Super-array to store compressed copies of the array arguments
real(kind=real_cvprec), allocatable :: args_cmpr(:,:)

! Super-array addresses
integer, parameter :: i_ref_temp   = 1
integer, parameter :: i_q_cond     = 2
integer, parameter :: i_rho_dry    = 3
integer, parameter :: i_rho_wet    = 4
integer, parameter :: i_dt_over_lz = 5
integer, parameter :: i_q_loc_cond = 6
integer, parameter :: i_n_cond     = 7
integer, parameter :: i_r_cond     = 8
integer, parameter :: i_wf_cond    = 9
integer, parameter :: i_kq_cond    = 10
integer, parameter :: i_kt_cond    = 11
integer, parameter :: n_args = 11

! Work variables for indexing null points for zero-reset
integer :: nc_null
integer, allocatable :: index_null(:)

! Loop counters
integer :: ic, ic2


if ( real(nc,real_cvprec) > cmpr_thresh * real(n_points,real_cvprec) ) then
  ! If majority of points have non-zero mixing ratio...

  ! Just call the main routine on the full arrays
  call calc_cond_properties( n_points,                                         &
                             cond_params, ref_temp, q_cond,                    &
                             rho_dry, rho_wet, dt_over_lz,                     &
                             q_loc_cond, n_cond,                               &
                             r_cond, wf_cond,                                  &
                             kq_cond, kt_cond )

  ! Reset outputs to zero at points where there is no condensate
  if ( nc < n_points ) then
    ! Store points where q_cond is zero
    allocate( index_null(n_points-nc) )
    nc_null = 0
    do ic = 1, n_points
      if ( .not. q_cond(ic) > zero ) then
        nc_null = nc_null + 1
        index_null(nc_null) = ic
      end if
    end do
    ! Reset outputs at these points
    do ic2 = 1, nc_null
      ic = index_null(ic2)
      q_loc_cond(ic) = zero
      n_cond(ic)     = zero
      r_cond(ic)     = zero
      wf_cond(ic)    = zero
      kq_cond(ic)    = zero
      kt_cond(ic)    = zero
    end do
    deallocate( index_null )
  end if


else  ! ( real(nc,real_cvprec) <= cmpr_thresh * real(n_points,real_cvprec) )
  ! If only a minority of points have non-zero mixing ratio...

  ! Do a compressed call...

  ! Allocate compression array
  allocate( args_cmpr ( nc, n_args  ) )

  ! Compress inputs
  do ic2 = 1, nc
    ic = index_ic(ic2)
    args_cmpr(ic2,i_ref_temp)   = ref_temp(ic)
    args_cmpr(ic2,i_q_cond)     = q_cond(ic)
    args_cmpr(ic2,i_rho_dry)    = rho_dry(ic)
    args_cmpr(ic2,i_rho_wet)    = rho_wet(ic)
    args_cmpr(ic2,i_dt_over_lz) = dt_over_lz(ic)
  end do

  ! Call the main routine on the compressed arrays
  call calc_cond_properties( nc,                                               &
         cond_params,                                                          &
         args_cmpr(:,i_ref_temp), args_cmpr(:,i_q_cond),                       &
         args_cmpr(:,i_rho_dry), args_cmpr(:,i_rho_wet),                       &
         args_cmpr(:,i_dt_over_lz),                                            &
         args_cmpr(:,i_q_loc_cond), args_cmpr(:,i_n_cond),                     &
         args_cmpr(:,i_r_cond), args_cmpr(:,i_wf_cond),                        &
         args_cmpr(:,i_kq_cond), args_cmpr(:,i_kt_cond) )

  ! Decompress outputs
  do ic2 = 1, nc
    ic = index_ic(ic2)
    q_loc_cond(ic) = args_cmpr(ic2,i_q_loc_cond)
    n_cond(ic)     = args_cmpr(ic2,i_n_cond)
    r_cond(ic)     = args_cmpr(ic2,i_r_cond)
    wf_cond(ic)    = args_cmpr(ic2,i_wf_cond)
    kq_cond(ic)    = args_cmpr(ic2,i_kq_cond)
    kt_cond(ic)    = args_cmpr(ic2,i_kt_cond)
  end do

  ! Deallocate
  deallocate( args_cmpr )

end if  ! ( real(nc,real_cvprec) <= cmpr_thresh * real(n_points,real_cvprec) )


return
end subroutine calc_cond_properties_cmpr

end module calc_cond_properties_mod
