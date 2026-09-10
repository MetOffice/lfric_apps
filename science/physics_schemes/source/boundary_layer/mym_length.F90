! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate mixing length in the MY model.
!           The square root of TKE, required by this subroutine and
!           elsewhere, is also returned

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  This code is based on the code provided by the authors who wrote
!  the following papers.
!   * Nakanishi, M. and H. Niino, 2009: Development of an improved
!      turbulence closure model for the atmospheric boundary layer.
!      J. Meteor. Soc. Japan, 87, 895-912.
!   * Nakanishi, M. and H. Niino, 2006: An improved Mellor-Yamada
!      Level-3 model: Its numerical stability and application to
!      a regional prediction of advection fog.
!      Boundary-Layer Meteor., 119, 397-407.
!   * Nakanishi, M. and H. Niino, 2004: An improved Mellor-Yamada
!      Level-3 model with condensation physics: Its design and
!      verification.
!      Boundary-Layer Meteor., 112, 1-31.
!   * Nakanishi, M., 2001: Improvement of the Mellor-Yamada
!      turbulence closure model based on large-eddy simulation data.
!      Boundary-Layer Meteor., 99, 349-378.
!   The web site publicising their code:
!    http://www.nda.ac.jp/~naka/MYNN/index.html

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_length_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_LENGTH_MOD'
contains

subroutine mym_length(                                                         &
      row_length, rows, halo_i, halo_j, bl_levels,                             &
      qke, z_uv, z_tq, dbdz, delta_smag, r_mosurf, fb_surf,                    &
      qkw, el)

use mym_const_mod, only: my_alpha4, one_third, elt_min, my_alpha1,             &
                         my_alpha2, my_alpha3
use mym_option_mod, only: tke_levels, my_z_limit_elb, l_3dtke, l_use_l_sq
use parkind1, only: jprb, jpim
use planet_constants_mod, only: vkman
use yomhook, only: lhook, dr_hook
use turb_diff_mod, only: mix_factor
implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   row_length,                                                                 &
                 ! Local number of points on a row
   rows,                                                                       &
                 ! Local number of rows in a theta field
   halo_i,                                                                     &
                 ! Size of halo in i direction.
   halo_j,                                                                     &
                 ! Size of halo in j direction.
   bl_levels
                 ! Max. no. of "boundary" levels

real(kind=r_bl), intent(in) ::                                          &
   qke(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                       &
                                                        bl_levels),            &
                 ! twice of TKE (denoted to q**2) on theta level K-1
   z_uv(row_length,rows,bl_levels+1),                                          &
                 ! Z_UV(*,K) is height of rho level k
   z_tq(row_length,rows,bl_levels),                                            &
                 ! Z_TQ(*,K) is height of theta level k.
   dbdz(row_length,rows,2:tke_levels),                                         &
                 ! Buoyancy gradient across layer
                 ! interface interpolated to theta levels.
                 ! (:,:,K) represents the value on theta level K-1
   delta_smag(row_length,rows),                                                &
                 ! IN delta_x used by Smagorinsky
   r_mosurf(row_length, rows),                                                 &
                 ! reciprocal of Monin-Obukhov Length
   fb_surf(row_length,rows)
                 ! Surface buoyancy flux over
                 ! density (m^2/s^3)

! Intent OUT Variables
real(kind=r_bl), intent(out) ::                                         &
   qkw(row_length, rows, tke_levels),                                          &
                 ! q=sqrt(qke) on theta level K-1
   el(row_length, rows, tke_levels)
                 ! mixing length on theta level K-1

! Local variables

integer ::                                                                     &
   i, j, k
                 ! Loop indexes
real(kind=r_bl) ::                                                      &
   qdz,                                                                        &
                 ! q times vertical grid space
   alp32,                                                                      &
                 ! combined constants (alpha3 / alpha2)
   rbv,                                                                        &
                 ! reciprocal of Brunt-Vaisala frequency
   elb,                                                                        &
                 ! mixing length related to buoyancy (L_B)
   els,                                                                        &
                 ! mixing length related to surface (L_S)
   ell,                                                                        &
                 ! additional mixing length for 3DTKE scheme (L_L)
   ellt,                                                                       &
                 ! temp mixing length for 3DTKE scheme
   zeta
                 ! non-dimensional length (height over MO length)

real(kind=r_bl) ::                                                      &
   elt(row_length, rows),                                                      &
                 ! mixing length related to vertical distribution
                 ! of TKE (L_T)
   vsc(row_length, rows)
                 ! work arrays

real(kind=r_bl), parameter ::                                           &
   zmax = 1.0,                                                                 &
                ! constant used in calculating els
   cns = 2.7
                ! constant used in calculating els

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_LENGTH'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do j = 1, rows
  do i = 1, row_length
    elt(i, j) = 0.0
    vsc(i, j) = 0.0
  end do
end do

do k = 1, tke_levels
  do j = 1, rows
    do i = 1, row_length
      qkw(i, j, k) = sqrt(max(qke(i, j, k), 1.0e-20))
    end do
  end do
end do

! vertical integration of qz and q
! Here, elt is still vertical integration of qz
! and vsc is that of q
do k = 2, tke_levels
  do j = 1, rows
    do i = 1, row_length
      qdz = qkw(i, j, k) * (z_uv(i, j, k) - z_uv(i, j, k - 1))
      elt(i, j) = elt(i, j) + qdz * z_tq(i, j, k - 1)
      vsc(i, j) = vsc(i, j) + qdz
    end do
  end do
end do

do j = 1, rows
  do i = 1, row_length
    elt(i, j) = max(my_alpha1 * elt(i, j) / (vsc(i, j) + 1.0e-10),             &
                    elt_min)
    vsc(i, j) = (elt(i, j) * max(fb_surf(i, j), 0.0)) ** one_third
  end do
end do

alp32 = my_alpha3 / my_alpha2
do k = 2, tke_levels
  do j = 1, rows
    do i = 1, row_length

      if (dbdz(i, j, k) > 0.0) then
        rbv = 1.0 / sqrt(dbdz(i, j, k))
        elb = my_alpha2 * qkw(i, j, k) * rbv                                   &
                  * (1.0 + alp32 * sqrt(vsc(i, j) * rbv / elt(i, j)))
      else
        elb = 1.0e10
      end if

      if (z_tq(i, j, k - 1) > my_z_limit_elb) then
        elb = min(elb, z_uv(i, j, k) - z_uv(i, j, k - 1))
      end if

      zeta = z_tq(i, j, k - 1) * r_mosurf(i, j)
      if (zeta > 0.0) then
        els = vkman * z_tq(i, j, k - 1)                                        &
                   / (1.0 + cns * min(zeta, zmax))
      else
        els = vkman * z_tq(i, j, k - 1)                                        &
                   * min((1.0 - my_alpha4 * zeta) ** 0.2, 2.0)
      end if
      if (l_3dtke) then
        ell = mix_factor * delta_smag(i,j)
        if (l_use_l_sq) then
          ellt = elb / ( elb / elt(i, j) + elb / els + 1.0)
          el(i, j, k) = sqrt(1.0/(1.0/(ellt*ellt) + 1.0/(ell*ell)))
        else
          el(i, j, k) = elb / ( elb / elt(i, j) + elb / els + elb / ell + 1.0)
        end if
      else
        el(i, j, k) = elb / ( elb / elt(i, j) + elb / els + 1.0)
      end if
    end do
  end do
end do

do j = 1, rows
  do i = 1, row_length
    el(i, j, 1) = el(i, j, 2)
  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_length
end module mym_length_mod
