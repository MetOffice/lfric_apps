! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate the mixing length in the first order closure
!           model

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module ddf_mix_length_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'DDF_MIX_LENGTH_MOD'
contains

subroutine ddf_mix_length(                                                     &
      row_length, rows, halo_i, halo_j, bl_levels,                             &
      z_uv, z_tq, dbdz, delta_smag, r_mosurf, fb_surf, h_pbl, e_trb,           &
      elm, coef_ce, ekw)

use mym_option_mod, only: tke_dlen,                                            &
            my_length, ddf_length, non_local_like_length,                      &
            l_tke_dlen_blackadar, tke_levels
use parkind1, only: jprb, jpim
use planet_constants_mod, only: vkman
use atm_fields_bounds_mod, only: tdims
use yomhook, only: lhook, dr_hook
use mym_length_mod, only: mym_length
implicit none

! Intent In Variables
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
   z_uv(row_length,rows,bl_levels+1),                                          &
                  ! Z_UV(*,K) is height of u level k
   z_tq(row_length,rows,bl_levels),                                            &
                  ! IN Z_TQ(*,K) is height of theta level k.
                  ! Cloud ice (kg per kg air)
   dbdz(row_length,rows,tke_levels),                                           &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   delta_smag(row_length,rows),                                                &
                  ! IN delta_x used by Smagorinsky
   r_mosurf(row_length, rows),                                                 &
                  ! reciprocal of Monin-Obkhov length
   fb_surf(row_length,rows),                                                   &
                  ! Surface flux buoyancy over density (m^2/s^3)
   h_pbl(row_length, rows),                                                    &
                  ! height of PBL determined by vertical profile
                  ! of SL
   e_trb(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end, bl_levels)
                  ! TKE defined on theta levels K-1

real(kind=r_bl), intent(out) ::                                         &
   elm(row_length, rows, tke_levels),                                          &
                  ! mixing length
   coef_ce(row_length, rows, tke_levels),                                      &
                  ! coefficient appeared in a dissipation term
   ekw(row_length, rows, tke_levels)
                  ! SQRT(e_trb)

! Local variables
integer :: i, j, k
                 ! loop counter

real(kind=r_bl) ::                                                      &
   rbv,                                                                        &
                  ! reciprocal of Brunt-Vaisala frequency
   elb,                                                                        &
                  ! mixing length driven by buoyancy
   els,                                                                        &
                  ! mixing length driven by surface
   delta_z
                  ! vertical grid spacing

real(kind=r_bl) ::                                                      &
   qke(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                       &
                                                  bl_levels),                  &
                  ! twice of TKE (denoted to q**2) on theta level K-1
   qkw(row_length, rows, tke_levels)
                  ! q=sqrt(qke) on theta level K-1

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DDF_MIX_LENGTH'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 2, tke_levels
  do j = 1, rows
    do i = 1, row_length
      ekw(i, j, k) = sqrt(max(e_trb(i, j, k), 1.0e-20))
    end do
  end do
end do

if (tke_dlen == my_length) then
  do k = 2, tke_levels
    do j = 1, rows
      do i = 1, row_length
        qke(i, j, k) = 2.0 * e_trb(i, j, k)
      end do
    end do
  end do
  call mym_length(                                                             &
        row_length, rows, halo_i, halo_j, bl_levels,                           &
        qke, z_uv, z_tq, dbdz, delta_smag, r_mosurf, fb_surf,                  &
        qkw, elm)
else if (tke_dlen == ddf_length                                                &
   .or. tke_dlen == non_local_like_length) then
  do k = 2, tke_levels
    do j = 1, rows
      do i = 1, row_length
        delta_z = z_uv(i, j, k) - z_uv(i, j, k - 1)
        if (dbdz(i, j, k) > 0.0) then
          rbv = 1.0 / sqrt(dbdz(i, j, k))
          elb = max(min(0.76 * ekw(i, j, k) * rbv,                             &
                        delta_z), 1.0e-10)
        else
          elb = delta_z
        end if
        elm(i, j, k) = elb
      end do
    end do
  end do

  if (tke_dlen == non_local_like_length) then
    do k = 2, tke_levels
      do j = 1, rows
        do i = 1, row_length
          if (z_tq(i, j, k - 1) < h_pbl(i, j) ) then
            elm(i, j, k) = 0.25 * 1.8 * h_pbl(i, j)                            &
                   * (1.0 - exp(                                               &
                             -4.0 * z_tq(i, j, k - 1)/h_pbl(i, j))             &
                      - 0.0003 * exp(                                          &
                            8.0 * z_tq(i, j, k - 1) / h_pbl(i, j)))
          end if
        end do
      end do
    end do
  end if  ! if tke_dlen == non_local_like_length

  if (l_tke_dlen_blackadar) then
    do k = 2, tke_levels
      do j = 1, rows
        do i = 1, row_length
          els = vkman * z_tq(i, j, k - 1)
          elm(i, j, k) = els / (1.0 + els / elm(i, j, k))
        end do
      end do
    end do
  end if
end if

! for diagnostics
do j = 1, rows
  do i = 1, row_length
    elm(i, j, 1) = elm(i, j, 2)
  end do
end do

if (tke_dlen == non_local_like_length) then
  do k = 2, tke_levels
    do j = 1, rows
      do i = 1, row_length
        coef_ce(i, j, k) = 0.41
      end do
    end do
  end do
else
  do k = 2, tke_levels
    do j = 1, rows
      do i = 1, row_length
        coef_ce(i, j, k) = 0.19 + 0.74 * elm(i, j, k)                          &
                              / (z_uv(i, j, k) - z_uv(i, j, k - 1))
      end do
    end do
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ddf_mix_length
end module ddf_mix_length_mod
