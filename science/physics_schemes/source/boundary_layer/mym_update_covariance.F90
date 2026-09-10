! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To integrate the covariances(tsq, qsq, cov) appeared
!           in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_update_covariance_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_UPDATE_COVARIANCE_MOD'
contains

subroutine mym_update_covariance(                                              &
! IN levels
      bl_levels,                                                               &
! IN fields
      z_uv, z_tq,                                                              &
      qkw, el, dfm, pdt_tsq, pdt_cov, pdt_res,                                 &
      pdq_qsq, pdq_cov, pdq_res, pdc_cov, pdc_tsq, pdc_qsq, pdc_res,           &
! INOUT fields
      tsq, qsq, cov)

use atm_fields_bounds_mod, only: tdims, tdims_s, pdims
use mym_const_mod, only: b2, coef_trbvar_diff
use mym_option_mod, only: l_my_extra_level, tke_levels
use timestep_mod, only: timestep
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
use mym_diff_matcoef_mod, only: mym_diff_matcoef
use mym_solve_simeq_mod, only: mym_solve_simeq
implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   bl_levels
                 ! Max. no. of "boundary" level

real(kind=r_bl), intent(in) ::                                          &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels+1),                                                          &
                 ! Z_UV(*,K) is height of u level k
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                 ! Z_TQ(*,K) is height of theta level k.
   qkw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                 ! sqrt(qke) = sqrt(2TKE)
   el(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! mixing length
   dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,            &
                                                  bl_levels),                  &
                 ! diffusion coefficients fot momentum
   pdt_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to tsq in the production term of tsq
   pdt_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to cov in the production term of tsq
   pdt_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a residual part in the production term of tsq
   pdq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to qsq in the production term of qsq
   pdq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to cov in the production term of qsq
   pdq_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a residual part in the production term of qsq
   pdc_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to cov in the production term of cov
   pdc_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to tsq in the production term of cov
   pdc_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to qsq in the production term of cov
   pdc_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels)
                 ! a residual part in the production term of cov

real(kind=r_bl), intent(in out) ::                                      &
   tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       bl_levels),                                                             &
                 ! Self covariance of liquid potential temperature
                 ! (thetal'**2) defined on theta levels K-1
   qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       bl_levels),                                                             &
                 ! Self covariance of total water
                 ! (qw'**2) defined on theta levels K-1
   cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       bl_levels)
                 ! Correlation between thetal and qw
                 ! (thetal'qw') defined on theta levels K-1

! Local Variables
integer ::                                                                     &
   i, j, k, k_start
                ! loop indexes, etc.
real(kind=r_bl) ::                                                      &
   elem
                ! work variables

real(kind=r_bl) ::                                                      &
   disp_coef
                ! coefficients of the prognostic variables in
                ! dissipation terms

real(kind=r_bl) ::                                                      &
   aa(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   bb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
   cc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                ! tri-diagonal matrix elements due to diffusion
   qq_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   qq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   qq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   aa_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   bb_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   cc_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   pp_tc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
   aa_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   bb_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   cc_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   pp_qc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
   aa_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   bb_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   cc_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          tke_levels),                                                         &
   pp_ct(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
   pp_cq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels)
                ! matrix elements (see the documents for details)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_UPDATE_COVARIANCE'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call mym_diff_matcoef(                                                         &
      bl_levels,coef_trbvar_diff, z_uv, z_tq, dfm, aa, bb, cc)

if (l_my_extra_level) then
  k_start = 1
else
  k_start = 2
end if

! set maxtrix elements
do k = k_start, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      disp_coef = 2.0 * qkw(i, j, k) / (b2 * el(i, j, k))
      elem = 1.0 - bb(i, j, k) * timestep                                      &
                             + timestep * disp_coef
      bb_tsq(i, j, k) = elem                                                   &
                          - 2.0 * pdt_tsq(i, j, k) * timestep
      bb_qsq(i, j, k) = elem                                                   &
                          - 2.0 * pdq_qsq(i, j, k) * timestep
      bb_cov(i, j, k) = elem                                                   &
                          - 2.0 * pdc_cov(i, j, k) * timestep
      qq_tsq(i, j, k) = tsq(i, j, k)                                           &
                        + 2.0 * pdt_res(i, j, k) * timestep
      qq_qsq(i, j, k) = qsq(i, j, k)                                           &
                        + 2.0 * pdq_res(i, j, k) * timestep
      qq_cov(i, j, k) = cov(i, j, k)                                           &
                        + 2.0 * pdc_res(i, j, k) * timestep

      elem = -aa(i, j, k) * timestep
      aa_tsq(i, j, k) = elem
      aa_qsq(i, j, k) = elem
      aa_cov(i, j, k) = elem

      elem = -cc(i, j, k) * timestep
      cc_tsq(i, j, k) = elem
      cc_qsq(i, j, k) = elem
      cc_cov(i, j, k) = elem

      pp_tc(i, j, k) = - 2.0 * pdt_cov(i, j, k) * timestep
      pp_qc(i, j, k) = - 2.0 * pdq_cov(i, j, k) * timestep
      pp_ct(i, j, k) = - 2.0 * pdc_tsq(i, j, k) * timestep
      pp_cq(i, j, k) = - 2.0 * pdc_qsq(i, j, k) * timestep
    end do
  end do
end do

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    aa_tsq(i, j, k_start) = 0.0
    aa_qsq(i, j, k_start) = 0.0
    aa_cov(i, j, k_start) = 0.0

    cc_tsq(i, j, tke_levels) = 0.0
    cc_qsq(i, j, tke_levels) = 0.0
    cc_cov(i, j, tke_levels) = 0.0
  end do
end do

if (.not. l_my_extra_level) then
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      bb_tsq(i, j, 1) = 1.0
      bb_qsq(i, j, 1) = 1.0
      bb_cov(i, j, 1) = 1.0
      qq_tsq(i, j, 1) = 0.0
      qq_qsq(i, j, 1) = 0.0
      qq_cov(i, j, 1) = 0.0

      aa_tsq(i, j, 1) = 0.0
      aa_qsq(i, j, 1) = 0.0
      aa_cov(i, j, 1) = 0.0

      cc_tsq(i, j, 1) = 0.0
      cc_qsq(i, j, 1) = 0.0
      cc_cov(i, j, 1) = 0.0

      pp_tc(i, j, 1) = 0.0
      pp_qc(i, j, 1) = 0.0
      pp_ct(i, j, 1) = 0.0
      pp_cq(i, j, 1) = 0.0
    end do
  end do
end if

! Solve the simultaneous equations for tsq, qsq and cov
call mym_solve_simeq(                                                          &
! IN levels
        bl_levels,                                                             &
! IN fields
        qq_tsq, qq_qsq, qq_cov, aa_tsq, bb_tsq, cc_tsq, pp_tc,                 &
        aa_qsq, bb_qsq, cc_qsq, pp_qc,aa_cov,bb_cov,cc_cov,pp_ct, pp_cq,       &
! OUT fields
        tsq, qsq, cov)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_update_covariance
end module mym_update_covariance_mod
