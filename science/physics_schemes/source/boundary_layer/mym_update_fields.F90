! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To integrate the prognostic variables appearing
!           in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_update_fields_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_UPDATE_FIELDS_MOD'
contains

subroutine mym_update_fields(bl_levels,coef,z_uv,z_tq,dfm,prod,disp_coef,field)

use atm_fields_bounds_mod, only: pdims, pdims_l, tdims, tdims_s
use mym_option_mod, only: l_my_extra_level, tke_levels
use timestep_mod, only: timestep
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
use mym_diff_matcoef_mod, only: mym_diff_matcoef
use mym_implic_mod, only: mym_implic
implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   bl_levels
                 ! Max. no. of "boundary" levels

real(kind=r_bl), intent(in) ::                                          &
   coef
                 ! factor for the diffusion coefficients to those for
                 ! momentum

real(kind=r_bl), intent(in) ::                                          &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
       bl_levels+1),                                                           &
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
       bl_levels),                                                             &
                 ! Z_TQ(*,K) is height of theta level k
   dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,            &
       bl_levels),                                                             &
                 ! diffusion coefficients for momentum
   prod(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! production term
   disp_coef(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tke_levels)
                 ! coefficients of dissipation term

! Intent INOUT Variables
real(kind=r_bl), intent(in out) ::                                      &
   field(pdims_l%i_start:pdims_l%i_end,pdims_l%j_start:pdims_l%j_end,          &
         bl_levels)
                 ! field to be integrated

! Local variables
integer ::                                                                     &
   i, j, k, k_start
                 ! Loop indexes

real(kind=r_bl) ::                                                      &
   aa(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),         &
   bb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),         &
   cc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),         &
   qq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels)
                ! coefficients of tri-diagonal equations

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_UPDATE_FIELDS'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the coefficients of tri-diagonal eqs. due to diffusion
call mym_diff_matcoef(bl_levels, coef, z_uv, z_tq, dfm, aa, bb, cc)

if (l_my_extra_level) then
  k_start = 1
else
  k_start = 2
end if

do k = k_start, tke_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      aa(i, j, k) = - aa(i, j, k) * timestep
      bb(i, j, k) = 1.0 - bb(i, j, k) * timestep                               &
                          + timestep * disp_coef(i, j, k)
      cc(i, j, k) = - cc(i, j, k) * timestep
      qq(i, j, k) = field(i, j, k) + timestep * prod(i, j, k)
    end do
  end do
end do

! Solve the tri-diagonal equations
call mym_implic(tke_levels, k_start, tke_levels, aa, bb, cc, qq)

do k = k_start, tke_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      field(i, j, k) = qq(i, j, k)
    end do
  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mym_update_fields
end module mym_update_fields_mod
