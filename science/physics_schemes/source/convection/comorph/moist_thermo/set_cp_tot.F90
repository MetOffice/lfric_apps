! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_cp_tot_mod

implicit none

contains

! Subroutine to compute the total heat capacity per unit
! dry-mass, as a function of the mixing ratios of all the
! water species

!----------------------------------------------------------------
! Version for compressed arrays with all the condensed water
! species in a super-array
!----------------------------------------------------------------
subroutine set_cp_tot( n_points, n_points_super,                               &
                       q_vap, q_cond_super, cp_tot )

use comorph_constants_mod, only: real_cvprec, n_cond_species,                  &
                     cp_dry, cp_vap, cond_params, zero

implicit none

! Number of points where we actually want to do something
integer, intent(in) :: n_points
! Number of points in the compressed super-array
! (to avoid having to repeatedly deallocate and reallocate when
!  re-using the same array, it maybe bigger than needed here).
integer, intent(in) :: n_points_super

! Water vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Super-array containing all condensed water mixing ratios
real(kind=real_cvprec), intent(in) :: q_cond_super                             &
                              ( n_points_super, n_cond_species )

! Parcel total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(out) :: cp_tot(n_points)

! Loop counters
integer :: i, i_cond

! Initialise to zero
do i = 1, n_points
  cp_tot(i) = zero
end do

! Loop over all condensed water species
do i_cond = 1, n_cond_species
  do i = 1, n_points
    ! Add contribution of each condensed water species onto
    ! the total heat capacity, using its mixing ratio
    ! and its specific heat capacity stored in cond_params
    cp_tot(i) = cp_tot(i) + q_cond_super(i,i_cond)                             &
                            * cond_params(i_cond)%pt % cp
  end do
end do

! Note: the larger contributions from dry air and vapour
! are added on last, because this leads to better precision
! in the summing of contributions from condensed water
! species above.

! Add on contributions from dry air and water vapour
do i = 1, n_points
  cp_tot(i) = cp_tot(i) + cp_dry + q_vap(i)*cp_vap
end do


return
end subroutine set_cp_tot


!----------------------------------------------------------------
! Version for 2-D slices of the full arrays
!----------------------------------------------------------------
subroutine set_cp_tot_2d( lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,                &
                          lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,                &
                          lb_s, ub_s, q_snow, lb_g, ub_g, q_graup,             &
                          cp_tot )

use comorph_constants_mod, only: real_hmprec,                                  &
                     l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,                &
                     cp_dry, cp_vap, cp_liq, cp_ice,                           &
                     nx_full, ny_full

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

! Mixing ratios of water species
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

! Parcel total heat capacity per unit dry-mass
real(kind=real_hmprec), intent(out) :: cp_tot(nx_full,ny_full)

! Copies of the heat capacities converted to host-model precision
real(kind=real_hmprec) :: cp_dry_p, cp_vap_p, cp_liq_p, cp_ice_p

integer :: i, j                   ! Loop counters

! Convert constants to same precision as the full 3-D fields
cp_dry_p = real(cp_dry,real_hmprec)
cp_vap_p = real(cp_vap,real_hmprec)
cp_liq_p = real(cp_liq,real_hmprec)
cp_ice_p = real(cp_ice,real_hmprec)

do j = 1, ny_full
  do i = 1, nx_full
    cp_tot(i,j) = q_cl(i,j)*cp_liq_p
  end do
end do

if ( l_cv_rain ) then
  do j = 1, ny_full
    do i = 1, nx_full
      cp_tot(i,j) = cp_tot(i,j) + q_rain(i,j)*cp_liq_p
    end do
  end do
end if

if ( l_cv_cf ) then
  do j = 1, ny_full
    do i = 1, nx_full
      cp_tot(i,j) = cp_tot(i,j) + q_cf(i,j)*cp_ice_p
    end do
  end do
end if

if ( l_cv_snow ) then
  do j = 1, ny_full
    do i = 1, nx_full
      cp_tot(i,j) = cp_tot(i,j) + q_snow(i,j)*cp_ice_p
    end do
  end do
end if

if ( l_cv_graup ) then
  do j = 1, ny_full
    do i = 1, nx_full
      cp_tot(i,j) = cp_tot(i,j) + q_graup(i,j)*cp_ice_p
    end do
  end do
end if

do j = 1, ny_full
  do i = 1, nx_full
    cp_tot(i,j) = cp_tot(i,j) + cp_dry_p + q_vap(i,j)*cp_vap_p
  end do
end do

return
end subroutine set_cp_tot_2d


end module set_cp_tot_mod
