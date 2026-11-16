! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To solve simultaneous equations for tsq, qsq and cov

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_solve_simeq_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_SOLVE_SIMEQ_MOD'
contains

subroutine mym_solve_simeq(                                                    &
! IN levels
      bl_levels,                                                               &
! IN fields
      qq_tsq, qq_qsq, qq_cov, aa_tsq, bb_tsq, cc_tsq, pp_tc,                   &
      aa_qsq, bb_qsq, cc_qsq, pp_qc,aa_cov,bb_cov, cc_cov, pp_ct, pp_cq,       &
! OUT fields
      tsq, qsq, cov)

use atm_fields_bounds_mod, only: tdims
use mym_option_mod, only: tke_levels, simeq_solver, bicgstab, gauss
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
use mym_solve_simeq_bcgstab_mod, only: mym_solve_simeq_bcgstab
use mym_solve_simeq_lud_mod, only: mym_solve_simeq_lud
use mym_solve_simeq_gauss_mod, only: mym_solve_simeq_gauss
implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   bl_levels
                 ! Max. no. of "boundary" level

real(kind=r_bl), intent(in) ::                                          &
   ! matrix elements (for meanings of each, see the document)
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
          tke_levels),                                                         &
   pp_cq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels)

real(kind=r_bl), intent(out) ::                                         &
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

! Local variables
integer ::                                                                     &
   i, j, k,                                                                    &
   endflag

real(kind=r_bl) ::                                                      &
   ! one-dimensional variables to secure continuous memory accesses
   qq_tsq_k(tke_levels),                                                       &
   qq_qsq_k(tke_levels),                                                       &
   qq_cov_k(tke_levels),                                                       &
   aa_tsq_k(tke_levels),                                                       &
   bb_tsq_k(tke_levels),                                                       &
   cc_tsq_k(tke_levels),                                                       &
   pp_tc_k(tke_levels),                                                        &
   aa_qsq_k(tke_levels),                                                       &
   bb_qsq_k(tke_levels),                                                       &
   cc_qsq_k(tke_levels),                                                       &
   pp_qc_k(tke_levels),                                                        &
   aa_cov_k(tke_levels),                                                       &
   bb_cov_k(tke_levels),                                                       &
   cc_cov_k(tke_levels),                                                       &
   pp_ct_k(tke_levels),                                                        &
   pp_cq_k(tke_levels),                                                        &
   tsq_k(tke_levels),                                                          &
   qsq_k(tke_levels),                                                          &
   cov_k(tke_levels)

! Parameters
integer, parameter ::                                                          &
   max_itr    = 500
             ! the maximum iteration number

real(kind=r_bl), parameter ::                                           &
   eps       = 1.0e-15
             ! convergence creteria

real(kind=r_bl), parameter ::                                           &
   tsq_scale = 1.0e0,                                                          &
   qsq_scale = 1.0e6,                                                          &
   cov_scale = 1.0e3,                                                          &
   r_tsq_scale = 1.0 / tsq_scale,                                              &
   r_qsq_scale = 1.0 / qsq_scale,                                              &
   r_cov_scale = 1.0 / cov_scale,                                              &
   tc_scale = tsq_scale * r_cov_scale,                                         &
   qc_scale = qsq_scale * r_cov_scale,                                         &
   ct_scale = cov_scale * r_tsq_scale,                                         &
   cq_scale = cov_scale * r_qsq_scale
             ! scaling factors for the matrix elements

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_SOLVE_SIMEQ'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end

    ! Copy to 1dim variables to secure continuous memory accesses
    do k = 1, tke_levels
      qq_tsq_k(k) = qq_tsq(i, j, k) * tsq_scale
      qq_qsq_k(k) = qq_qsq(i, j, k) * qsq_scale
      qq_cov_k(k) = qq_cov(i, j, k) * cov_scale
      aa_tsq_k(k) = aa_tsq(i, j, k)
      bb_tsq_k(k) = bb_tsq(i, j, k)
      cc_tsq_k(k) = cc_tsq(i, j, k)
      pp_tc_k(k)  = pp_tc(i, j, k)  * tc_scale
      aa_qsq_k(k) = aa_qsq(i, j, k)
      bb_qsq_k(k) = bb_qsq(i, j, k)
      cc_qsq_k(k) = cc_qsq(i, j, k)
      pp_qc_k(k)  = pp_qc(i, j, k)  * qc_scale
      aa_cov_k(k) = aa_cov(i, j, k)
      bb_cov_k(k) = bb_cov(i, j, k)
      cc_cov_k(k) = cc_cov(i, j, k)
      pp_ct_k(k)  = pp_ct(i, j, k)  * ct_scale
      pp_cq_k(k)  = pp_cq(i, j, k)  * cq_scale
    end do

    if (simeq_solver == bicgstab) then
      call  mym_solve_simeq_bcgstab(                                           &
              max_itr, eps,                                                    &
              qq_tsq_k, qq_qsq_k, qq_cov_k,                                    &
              aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                           &
              aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                           &
              aa_cov_k, bb_cov_k, cc_cov_k,                                    &
              pp_ct_k, pp_cq_k,                                                &
              tsq_k, qsq_k, cov_k, endflag)

      if (endflag < 0) then
        ! if failed to converge, solve eqs. by LU decomposition
        call mym_solve_simeq_lud(                                              &
              qq_tsq_k, qq_qsq_k, qq_cov_k,                                    &
              aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                           &
              aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                           &
              aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                  &
              tsq_k, qsq_k, cov_k)
      end if
    else if (simeq_solver == gauss) then
      call mym_solve_simeq_gauss(                                              &
              qq_tsq_k, qq_qsq_k, qq_cov_k,                                    &
              aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                           &
              aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                           &
              aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                  &
              tsq_k, qsq_k, cov_k)
    end if

    ! set the values into the original arrays.
    do k = 1, tke_levels
      tsq(i, j, k) = tsq_k(k) * r_tsq_scale
      qsq(i, j, k) = qsq_k(k) * r_qsq_scale
      cov(i, j, k) = cov_k(k) * r_cov_scale
    end do

  end do !loop i = tdims%i_start, tdims%i_end
end do   !loop j = tdims%j_start, tdims%j_end
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_solve_simeq
end module mym_solve_simeq_mod
