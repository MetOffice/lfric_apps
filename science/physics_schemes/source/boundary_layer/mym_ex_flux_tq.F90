! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To calculate heat and moisture fluxes in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_ex_flux_tq_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_EX_FLUX_TQ_MOD'
contains

subroutine mym_ex_flux_tq(                                                     &
      bl_levels,                                                               &
      tl, qw, rhokh, rhogamt, rhogamq, rdz,                                    &
      ftl, fqw)

use atm_fields_bounds_mod, only: tdims, pdims
use model_domain_mod,      only: model_type, mt_single_column
use planet_constants_mod,  only: cp, grcp

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! INTENT IN Variables
integer, intent(in) ::                                                         &
   bl_levels
                 ! Max. no. of "boundary" levels

real(kind=r_bl), intent(in) ::                                          &
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                   ! Liquid/frozen water temperture (K)
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),         &
                   ! Total water content (kg/kg)
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels),                                                           &
                   ! Exchange coeffs for scalars
                   ! between K and K-1 on theta levels.
                   ! i.e. the coeffs are defined on rho levels
   rhogamt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels),                                                       &
                   ! Counter gradient term for FTL on rho levels
   rhogamq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels),                                                       &
                   ! Counter gradient Term for FQW on
   rdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)
                   ! RDZ(,1) is the reciprocal
                   ! height of level 1, i.e. of the
                   ! middle of layer 1.  For K > 1,
                   ! RDZ(,K) is the reciprocal of the
                   ! vertical distance from level
                   ! K-1 to level K.

! INTENT OUT Variables
real(kind=r_bl), intent(in out) ::                                      &
   ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),         &
                   ! FTL(,K) contains net turb
                   ! sensible heat flux into layer K
                   ! from below; so FTL(,1) is the
                   ! surface sensible heat, H. (W/m2)
                   ! defined on rho levels
   fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, bl_levels)
                   ! Moisture flux between layers
                   ! (kg per square metre per sec).
                   ! FQW(,1) is total water flux
                   ! from surface, 'E'.
                   ! defined on rho levels

character(len=*), parameter ::  RoutineName = 'MYM_EX_FLUX_TQ'

! LOCAL VARIABLES.

integer ::                                                                     &
   i, j, k

real(kind=r_bl) ::                                                      &
   grad_ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                   ! Gradient part of FTL
                   ! K*dth/dz
   grad_fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                   ! Gradient part of FQW
                   ! K*dq/dz
   count_grad_ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  bl_levels),                                                  &
                   ! Counter gradient part of FTL
   count_grad_fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  bl_levels)
                   ! Counter gradient part of FQW

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 1, bl_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      grad_ftl(i,j,k)=0.0
      grad_fqw(i,j,k)=0.0
      count_grad_ftl(i,j,k)=0.0
      count_grad_fqw(i,j,k)=0.0
    end do
  end do
end do

do k = 2, bl_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      grad_ftl(i,j,k)= - rhokh(i,j,k) *                                        &
                ( ( ( tl(i,j,k) - tl(i,j,k-1) ) * rdz(i,j,k) )                 &
                                                        + grcp )
      grad_fqw(i,j,k)= - rhokh(i,j,k) *                                        &
                    ( qw(i,j,k) - qw(i,j,k-1) ) * rdz(i,j,k)
      count_grad_ftl(i,j,k) = -rhogamt(i,j,k)
      count_grad_fqw(i,j,k) = -rhogamq(i,j,k)
      ftl(i,j,k) = grad_ftl(i,j,k) + count_grad_ftl(i,j,k)
      fqw(i,j,k) = grad_fqw(i,j,k) + count_grad_fqw(i,j,k)
    end do
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mym_ex_flux_tq
end module mym_ex_flux_tq_mod
