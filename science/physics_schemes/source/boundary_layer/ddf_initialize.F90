! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To set the initial TKE in the first order closure model

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module ddf_initialize_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'DDF_INITIALIZE_MOD'
contains

subroutine ddf_initialize(                                                     &
      bl_levels,                                                               &
      z_uv, z_tq, dbdz, dvdzm, delta_smag, r_mosurf, fb_surf, u_s, h_pbl,      &
      e_trb)

use atm_fields_bounds_mod, only: tdims, pdims, tdims_s
use mym_const_mod, only: e_trb_max
use mym_option_mod, only: tke_levels, l_my_extra_level,                        &
                          my_z_extra_fact, my_lowest_pd_surf,                  &
                          tke_cm_mx, tke_cm_fa
use parkind1, only: jprb, jpim
use planet_constants_mod, only: vkman
use yomhook, only: lhook, dr_hook
use ddf_mix_length_mod, only: ddf_mix_length
use mym_calcphi_mod, only: mym_calcphi
use mym_diff_matcoef_mod, only: mym_diff_matcoef
use mym_implic_mod, only: mym_implic
implicit none

! Intent In Variables
integer, intent(in) ::                                                         &
   bl_levels
                  ! Max. no. of "boundary" levels

real(kind=r_bl), intent(in) ::                                          &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels+1),                                                          &
                  ! Z_UV(*,K) is height of u level k
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                  ! IN Z_TQ(*,K) is height of theta level k.
                  ! Cloud ice (kg per kg air)
   dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                  ! Modulus of wind shear at theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   delta_smag(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                  ! IN delta_x used by Smagorinsky
   r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                  ! reciprocal of Monin-Obkhov length
   fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
                  ! Surface flux buoyancy over density (m^2/s^3)
   u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                  ! Surface friction velocity
   h_pbl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                  ! height of PBL determined by vertical profile
                  ! of SL

real(kind=r_bl), intent(out) ::                                         &
   e_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
                                                  bl_levels)
                  ! TKE defined on theta levels K-1

! Local variables
integer ::                                                                     &
   i, j, k, ll,                                                                &
   itr_ini

real(kind=r_bl) ::                                                      &
   r_pr,                                                                       &
   elq,                                                                        &
   sm,                                                                         &
   sh,                                                                         &
   gm,                                                                         &
   gh

real(kind=r_bl) ::                                                      &
   ekw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
   elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
   coef_cm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
   coef_ce(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
   pdk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
   dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,            &
       tke_levels),                                                            &
   aa(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   bb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   cc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   pdk0(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
   pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
   phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

real(kind=r_bl), parameter ::                                           &
   pr = 0.7,                                                                   &
                  ! Prandtl number
                  ! only in the initialization,
                  ! constant prandtl number is assumed.
   diff_fact = 2.0
                  ! factor of a diffusion coef of E_TRB to that of
                  ! momentum

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DDF_INITIALIZE'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

r_pr = 1.0 / pr

if (my_lowest_pd_surf == 0) then
  l_my_extra_level = .false.
  my_z_extra_fact = 1.0
end if

! initial guess for e_trb, assuming neutral layer
! and set some parameters
do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if (z_tq(i, j, k - 1) < h_pbl(i, j)) then
        coef_cm(i, j, k) = tke_cm_mx
      else
        coef_cm(i, j, k) = tke_cm_fa
      end if
      sm = coef_cm(i, j, k)
      sh = coef_cm(i, j, k) * r_pr
      gm = dvdzm(i, j, k) ** 2
      gh = -dbdz(i, j, k)
      pdk(i, j, k) = sm * gm + sh * gh
      if (pdk(i, j, k) <= 0.0) then
        pdk(i, j, k) = 0.0
        e_trb(i, j, k) = 0.0
      else
        e_trb(i, j, k) = 1.0e-5
      end if
    end do
  end do
end do

if (my_lowest_pd_surf > 0) then
  call mym_calcphi(                                                            &
        bl_levels, z_tq, r_mosurf, pmz, phh)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      pdk0(i, j) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)                            &
         / (vkman * z_tq(i, j, 1))
    end do
  end do
end if  ! IF MY_lowest_pd_surf

itr_ini = tke_levels + 1

do ll = 1, itr_ini
  call ddf_mix_length(                                                         &
    tdims%i_end, tdims%j_end, 0, 0, bl_levels,                                 &
    z_uv, z_tq, dbdz, delta_smag, r_mosurf, fb_surf, h_pbl, e_trb,             &
    elm, coef_ce, ekw)

  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if (e_trb(i, j, k) <= 0.0) then
          ekw(i, j, k) = 0.0
        end if
        dfm(i, j, k) = coef_cm(i, j, k) * ekw(i, j, k) * elm(i, j, k)
      end do
    end do
  end do

  call mym_diff_matcoef(                                                       &
        bl_levels, diff_fact, z_uv, z_tq, dfm, aa, bb, cc)

  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if (bb(i, j, k) == 0.0) then
          aa(i, j, k) = 0.0
          bb(i, j, k) = 1.0
          cc(i, j, k) = 0.0
          e_trb(i, j, k) = 0.0
        else
          elq = ekw(i, j, k) * elm(i, j, k)
          aa(i, j, k) = - aa(i, j, k)
          bb(i, j, k) = - bb(i, j, k)                                          &
                  + ekw(i, j, k) * coef_ce(i, j, k)                            &
                                     / max(elm(i, j, k), 1.0e-20)
          bb(i, j, k) = sign(max(abs(bb(i, j, k)), 1.0e-20_r_bl),       &
                                    bb(i, j, k))

          cc(i, j, k) = - cc(i, j, k)
          e_trb(i, j, k) = elq * pdk(i, j, k)
        end if
      end do
    end do
  end do

  if (my_lowest_pd_surf > 0) then
    k = 2
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if (bb(i, j, k) /= 0.0 .and. pdk(i, j, k) > 0.0) then
          e_trb(i, j, k) = pdk0(i, j)
        end if
      end do
    end do
  end if   ! IF MY_lowest_pd_surf > 0

  call mym_implic(                                                             &
                  tke_levels, 2, tke_levels, aa, bb, cc, e_trb)
end do  ! DO ll = 1, itr_ini

do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      e_trb(i, j, k) = min(                                                    &
                            max(e_trb(i, j, k), 1.0e-20),                      &
                                e_trb_max)
    end do
  end do
end do

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    e_trb(i, j, 1) = 0.0
  end do
end do

do k = tke_levels + 1, bl_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      e_trb(i, j, k) = 0.0
    end do
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ddf_initialize
end module ddf_initialize_mod
