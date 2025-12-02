! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_q_tot_mod

implicit none

contains


!----------------------------------------------------------------
! Subroutines to calculate total-water mixing-ratio
! using the mixing ratios of all the water species.
!----------------------------------------------------------------
! Note that the dry-mass to wet-mass conversion factor is
! 1 + q_tot

! Version for compressed arrays, with the condensed water
! species in a super-array
subroutine calc_q_tot( n_points, n_points_super,                               &
                       q_vap, q_cond_super, q_tot )

use comorph_constants_mod, only: real_cvprec, n_cond_species, zero

implicit none

! Number of points where we actually want to do something
integer, intent(in) :: n_points
! Number of points in the compressed super-array
! (to avoid having to repeatedly deallocate and reallocate when
!  re-using the same array, it maybe bigger than needed here).
integer, intent(in) :: n_points_super

! Water vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Super-array containing all condensed water species
real(kind=real_cvprec), intent(in) :: q_cond_super                             &
                              ( n_points_super, n_cond_species )

! Output total-water mixing-ratio
real(kind=real_cvprec), intent(out) :: q_tot(n_points)

! Loop counters
integer :: i, i_cond

! Initialise to zero
do i = 1, n_points
  q_tot(i) = zero
end do

! Note: since q_vap is usually much larger than the condensed
! water fields, greater precision is achieved by adding up the
! condensed water fields first, then adding on the larger q_vap

! Add on condensed water mixing ratios
do i_cond = 1, n_cond_species
  do i = 1, n_points
    q_tot(i) = q_tot(i) + q_cond_super(i,i_cond)
  end do
end do

! Add on water vapour mixing ratio
do i = 1, n_points
  q_tot(i) = q_tot(i) +  q_vap(i)
end do

return
end subroutine calc_q_tot


!----------------------------------------------------------------
! Version for full 3-D arrays
!----------------------------------------------------------------
subroutine calc_q_tot_3d( lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,                &
                          lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,                &
                          lb_s, ub_s, q_snow, lb_g, ub_g, q_graup,             &
                          q_tot )

use comorph_constants_mod, only: real_hmprec,                                  &
                     l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,                &
                     nx_full, ny_full, k_bot_conv, k_top_conv

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

! Water vapour mixing ratio
integer, intent(in) :: lb_v(3), ub_v(3)
real(kind=real_hmprec), intent(in) :: q_vap                                    &
                         ( lb_v(1):ub_v(1), lb_v(2):ub_v(2), lb_v(3):ub_v(3) )

! Condensed water species mixing ratios
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

! Output dry-mass to wet-mass conversion factor
real(kind=real_hmprec), intent(out) :: q_tot                                   &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Loop counters
integer :: i, j, k

!$OMP PARALLEL DEFAULT(none) private( i, j, k )                                &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv,                        &
!$OMP         q_tot, q_cl, q_rain, q_cf, q_snow, q_graup, q_vap )

! Initialise using liquid cloud
!$OMP do SCHEDULE(STATIC)
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      q_tot(i,j,k) = q_cl(i,j,k)
    end do
  end do
end do
!$OMP end do

! Add on mixing ratios of optional hydrometeor species...

if ( l_cv_rain ) then
!$OMP do SCHEDULE(STATIC)
  do k = k_bot_conv, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        q_tot(i,j,k) = q_tot(i,j,k) + q_rain(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if

if ( l_cv_cf ) then
!$OMP do SCHEDULE(STATIC)
  do k = k_bot_conv, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        q_tot(i,j,k) = q_tot(i,j,k) + q_cf(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if

if ( l_cv_snow ) then
!$OMP do SCHEDULE(STATIC)
  do k = k_bot_conv, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        q_tot(i,j,k) = q_tot(i,j,k) + q_snow(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if

if ( l_cv_graup ) then
!$OMP do SCHEDULE(STATIC)
  do k = k_bot_conv, k_top_conv
    do j = 1, ny_full
      do i = 1, nx_full
        q_tot(i,j,k) = q_tot(i,j,k) + q_graup(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if

! Add on water-vapour
!$OMP do SCHEDULE(STATIC)
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      q_tot(i,j,k) = q_tot(i,j,k) + q_vap(i,j,k)
    end do
  end do
end do
!$OMP end do

!$OMP end PARALLEL

return
end subroutine calc_q_tot_3d


end module calc_q_tot_mod
