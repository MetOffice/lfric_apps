! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: To calculate gradient functions at the surface used in
!           evaluating the production terms of the prognostic
!           variables at the lowest layer.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_calcphi_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_CALCPHI_MOD'
contains

subroutine mym_calcphi(bl_levels, z_tq, r_mosurf, pmz, phh)

use atm_fields_bounds_mod, only: tdims
use mym_const_mod, only: two_thirds, pr
use mym_option_mod, only:                                                      &
      businger, bh1991, my_lowest_pd_surf,                                     &
      l_my_extra_level, my_z_extra_fact
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   bl_levels
                 ! number of boundary layer levels

real(kind=r_bl), intent(in) ::                                          &
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                 ! Z_TQ(*,K) is height of theta
                 !    level k.
   r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! reciprocal of Monin-Obkhov length

! Intent OUT Variables
real(kind=r_bl), intent(out) ::                                         &
   pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                 ! gradient function for momentum
                 ! at surface minus non-dimensional height
   phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! gradient function for scalars
                 ! at surface

! Local variables
integer ::                                                                     &
   i, j
                 ! Loop indexes

real(kind=r_bl) ::                                                      &
   zeta,                                                                       &
                 ! non-dimensional height
   tmp
                 ! work variable

real(kind=r_bl) ::                                                      &
   z_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! height of the lowest layer

real(kind=r_bl), parameter ::                                           &
                 ! coefficients appeared
                 !               in Beljaars and Holtslag(1991)
   bel_a = 1.0,                                                                &
   bel_b = 2.0 / 3.0,                                                          &
   bel_c = 5.0,                                                                &
   bel_d = 0.35

real(kind=r_bl), parameter ::                                           &
   my_zeta_max = 2.0
                ! upper limit for zeta

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_CALCPHI'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_my_extra_level) then
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      z_1(i, j) = z_tq(i, j, 1) * my_z_extra_fact
    end do
  end do
else
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      z_1(i, j) = z_tq(i, j, 1)
    end do
  end do
end if

if (my_lowest_pd_surf == businger) then
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      zeta = min(z_1(i, j) * r_mosurf(i, j), my_zeta_max)
      if (zeta >= 0.0) then
        pmz(i, j) = 1.0 + 4.7 * zeta
        phh(i, j) = pr + 4.7 * zeta
      else
        pmz(i, j) = 1.0 / sqrt(sqrt(1.0 - 15.0 * zeta))
        phh(i, j) = pr / sqrt(1.0 - 9.0 * zeta)
      end if
      pmz(i, j) = pmz(i, j) - zeta
    end do
  end do
else if (my_lowest_pd_surf == bh1991) then
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      zeta = min(z_1(i, j) * r_mosurf(i, j), my_zeta_max)
      if (zeta >= 0) then
        tmp = bel_b * exp(-bel_d * zeta)                                       &
             * (bel_d * zeta - bel_c - 1.0)
        pmz(i, j) = 1.0 - zeta * (tmp - bel_a)
        phh(i, j) = 1.0 - zeta * (tmp -                                        &
             sqrt(1.0 + two_thirds * bel_a * zeta))
      else
        tmp = sqrt(1.0 - 16.0 * zeta)
        pmz(i, j) = 1.0 / sqrt(tmp)
        phh(i, j) = 1.0 / tmp
      end if
      pmz(i, j) = pmz(i, j) - zeta
    end do
  end do
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_calcphi
end module mym_calcphi_mod
