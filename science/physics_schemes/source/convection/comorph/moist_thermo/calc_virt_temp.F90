! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_virt_temp_mod

implicit none

contains


!----------------------------------------------------------------
! Subroutine to calculate virtual temperature
!----------------------------------------------------------------

! Version for compressed arrays, with the condensed water
! species in a super-array
subroutine calc_virt_temp( n_points, n_points_super,                           &
                           temperature, q_vap, q_cond_super,                   &
                           virt_temp )

use comorph_constants_mod, only: real_cvprec, n_cond_species, one
use calc_virt_temp_dry_mod, only: calc_virt_temp_dry
use calc_q_tot_mod, only: calc_q_tot

implicit none

! Number of points where we actually want to do something
integer, intent(in) :: n_points
! Number of points in the compressed super-array
! (to avoid having to repeatedly deallocate and reallocate when
!  re-using the same array, it maybe bigger than needed here).
integer, intent(in) :: n_points_super

! Temperature
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Water vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Super-array containing all condensed water species
real(kind=real_cvprec), intent(in) :: q_cond_super                             &
                              ( n_points_super, n_cond_species )

! Output virtual temperature
real(kind=real_cvprec), intent(out) :: virt_temp(n_points)

! Work variable: total-water mixing-ratio
real(kind=real_cvprec) :: q_tot(n_points)

! Loop counters
integer :: i

! In terms of mixing ratios qv, qc, virtual temperature is:
! Tv = T ( 1 + Rv/Rd qv ) / ( 1 + qv + qc )

! Calculate the numerator T ( 1 + Rv/Rd qv ),
! which is the volume of the parcel per unit dry-mass
! * pressure/R_dry, = p/(R_dry rho_dry)
call calc_virt_temp_dry( n_points, temperature, q_vap,                         &
                         virt_temp )

! Calculate the total-water qv + qc
call calc_q_tot( n_points, n_points_super,                                     &
                 q_vap, q_cond_super, q_tot )

! Complete the formula for virtual temperature
do i = 1, n_points
  virt_temp(i) = virt_temp(i) / ( one + q_tot(i) )
end do

return
end subroutine calc_virt_temp


!----------------------------------------------------------------
! Version for full 3-D arrays
!----------------------------------------------------------------
subroutine calc_virt_temp_3d( lb_T, ub_T, temperature,                         &
                              lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,            &
                              lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,            &
                              lb_s, ub_s, q_snow, lb_g, ub_g, q_graup,         &
                              virt_temp )

use comorph_constants_mod, only: real_hmprec,                                  &
                     nx_full, ny_full, k_bot_conv, k_top_conv

use calc_virt_temp_dry_mod, only: calc_virt_temp_dry_3d
use calc_q_tot_mod, only: calc_q_tot_3d

implicit none

! Note: the array arguments below may have halos, and in here
! we don't care about the halos.
! We need to pass in the bounds of each array to use in
! the declarations, to ensure the do loops pick out the same
! indices of the arrays regardless of whether or not they
! have halos.  The lower and upper bounds are passed in through the
! argument list in the lb_* and ub_* integer arrays.  Those storing the
! bounds for 3D arrays must have 3 elements; one for each dimension
! of the array.

! Temperature
integer, intent(in) :: lb_T(3), ub_T(3)
real(kind=real_hmprec), intent(in) :: temperature                              &
                         ( lb_T(1):lb_T(1), lb_T(2):lb_T(2), lb_T(3):lb_T(3) )

! Water vapour mixing ratio
integer, intent(in) :: lb_v(3), ub_v(3)
real(kind=real_hmprec), intent(in) :: q_vap                                    &
                         ( lb_v(1):ub_v(1), lb_v(2):ub_v(2), lb_v(3):ub_v(3) )

! Hydrometeor mixing ratios
integer, intent(in) :: lb_l(3), ub_l(3)
real(kind=real_hmprec), intent(in) :: q_cl                                     &
                         ( lb_l(1):ub_l(1), lb_l(2):ub_l(2), lb_l(3):ub_l(3) )
integer, intent(in) :: lb_r(3), ub_r(3)
real(kind=real_hmprec), intent(in) :: q_rain                                   &
                         ( lb_r(1):ub_r(1), lb_r(2):ub_r(2), lb_r(3):ub_r(3) )
integer, intent(in) :: lb_f(3), ub_f(3)
real(kind=real_hmprec), intent(in) :: q_cf                                     &
                         ( lb_f(1):ub_f(1), lb_f(2):ub_f(2), lb_f(3):ub_f(3) )
integer, intent(in) :: lb_s(3), ub_s(3)
real(kind=real_hmprec), intent(in) :: q_snow                                   &
                         ( lb_s(1):ub_s(1), lb_s(2):ub_s(2), lb_s(3):ub_s(3) )
integer, intent(in) :: lb_g(3), ub_g(3)
real(kind=real_hmprec), intent(in) :: q_graup                                  &
                         ( lb_g(1):ub_g(1), lb_g(2):ub_g(2), lb_g(3):ub_g(3) )

! Output virtual temperature
real(kind=real_hmprec), intent(out) :: virt_temp                               &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Work variable: total-water mixing-ratio
real(kind=real_hmprec) :: q_tot                                                &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

real(kind=real_hmprec), parameter :: one = 1.0_real_hmprec

! Loop counters
integer :: i, j, k

! In terms of mixing ratios q, qc, virtual temperature is:
! Tv = T ( 1 + Rv/Rd qv ) / ( 1 + qv + qc )

! Calculate the numerator T ( 1 + Rv/Rd qv ),
! which is the volume of the parcel per unit dry-mass
! * pressure/R_dry, = p/(R_dry rho_dry)
call calc_virt_temp_dry_3d( lb_T, ub_T, temperature, lb_v, ub_v, q_vap,        &
                            virt_temp )

! Calculate the total-water qv + qc
call calc_q_tot_3d( lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,                      &
                    lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,                      &
                    lb_s, ub_s, q_snow, lb_g, ub_g, q_graup,                   &
                    q_tot )

! Complete the formula for virtual temperature
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC) private( i, j, k )            &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv, virt_temp, q_tot )
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      virt_temp(i,j,k) = virt_temp(i,j,k) / (one + q_tot(i,j,k))
    end do
  end do
end do
!$OMP end PARALLEL do

return
end subroutine calc_virt_temp_3d


end module calc_virt_temp_mod
