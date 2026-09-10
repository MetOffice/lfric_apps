! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To solve simultaneous equations by bi-conjugate gradient
!           stabilized method (BCGSTAB)
!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_solve_simeq_bcgstab_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter,                                                   &
                  private :: ModuleName = 'MYM_SOLVE_SIMEQ_BCGSTAB_MOD'
contains

subroutine mym_solve_simeq_bcgstab(                                            &
   max_itr, eps,                                                               &
   qq_tsq_k, qq_qsq_k, qq_cov_k,                                               &
   aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                      &
   aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                      &
   aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                             &
   tsq_k, qsq_k, cov_k, endflag)

use mym_option_mod, only: tke_levels
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use mym_simeq_ilud2_decmp_mod, only: mym_simeq_ilud2_decmp
use mym_simeq_matrix_prod_mod, only: mym_simeq_matrix_prod
use mym_solve_simeq_ilud2_mod, only: mym_solve_simeq_ilud2
implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   max_itr
             ! the maximum number of iterations

real(kind=r_bl), intent(in) ::                                          &
   eps
             ! convergence condition

real(kind=r_bl), intent(in) ::                                          &
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
   pp_cq_k(tke_levels)
           ! matrix elements

real(kind=r_bl), intent(out) ::                                         &
   tsq_k(tke_levels),                                                          &
   qsq_k(tke_levels),                                                          &
   cov_k(tke_levels)
           ! solved tsq, qsq and cov

integer, intent(out) ::                                                        &
   endflag
           ! to indicate if converged
           ! positive means proper solution is obtains.
           ! 0: converged
           ! 1: obtained an exact solution (residual = 0)
           ! -1: max_itr iterations were done, but not converged
           ! -2: solution in iterations becomes unexpectedly large,
           !     so gave up

! Local variables
integer ::                                                                     &
   k, m,                                                                       &
           ! loop indexes
   nitr
           ! a number of iterations

real(kind=r_bl) ::                                                      &
   norm,                                                                       &
           ! residual norm
   r_qq_norm,                                                                  &
           ! reciprocal of the inirial residual norm
   err,                                                                        &
           ! norm * r_qq_norm
   bet,                                                                        &
           ! beta
   alp_num,                                                                    &
           ! numerator of alpha
   alp_den,                                                                    &
           ! denominator of alpha
   alp,                                                                        &
           ! alpha
   omg_num,                                                                    &
           ! numerator of omega
   omg_den,                                                                    &
           ! denominator of omega
   omg,                                                                        &
           ! omega
   max_val
           ! maximum value of solutions

real(kind=r_bl) ::                                                      &
   rvec_tsq(tke_levels),                                                       &
   rvec_qsq(tke_levels),                                                       &
   rvec_cov(tke_levels),                                                       &
   r0vec_tsq(tke_levels),                                                      &
   r0vec_qsq(tke_levels),                                                      &
   r0vec_cov(tke_levels),                                                      &
   pvec_tsq(tke_levels),                                                       &
   pvec_qsq(tke_levels),                                                       &
   pvec_cov(tke_levels),                                                       &
   ppvec_tsq(tke_levels),                                                      &
   ppvec_qsq(tke_levels),                                                      &
   ppvec_cov(tke_levels),                                                      &
   vvec_tsq(tke_levels),                                                       &
   vvec_qsq(tke_levels),                                                       &
   vvec_cov(tke_levels),                                                       &
   svec_tsq(tke_levels),                                                       &
   svec_qsq(tke_levels),                                                       &
   svec_cov(tke_levels),                                                       &
   ssvec_tsq(tke_levels),                                                      &
   ssvec_qsq(tke_levels),                                                      &
   ssvec_cov(tke_levels),                                                      &
   tvec_tsq(tke_levels),                                                       &
   tvec_qsq(tke_levels),                                                       &
   tvec_cov(tke_levels),                                                       &
           ! Vectors used in the BCG algorithm.
           ! See the document
   aap_tsq_k(tke_levels),                                                      &
   r_bbp_tsq_k(tke_levels),                                                    &
   ccp_tsq_k(tke_levels),                                                      &
   aap_qsq_k(tke_levels),                                                      &
   r_bbp_qsq_k(tke_levels),                                                    &
   ccp_qsq_k(tke_levels),                                                      &
   aap_cov_k(tke_levels),                                                      &
   r_bbp_cov_k(tke_levels),                                                    &
   ccp_cov_k(tke_levels),                                                      &
   ppp_tc_k(tke_levels, 0:2),                                                  &
   ppp_qc_k(tke_levels, 0:2),                                                  &
   ppp_ct_k(tke_levels, 0:2),                                                  &
   ppp_cq_k(tke_levels, 0:2)
          ! elements of ILU(2)
          ! the second dimension corresponds to the fill-in level
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_SOLVE_SIMEQ_BCGSTAB'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call mym_simeq_ilud2_decmp(                                                    &
   aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                      &
   aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                      &
   aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                             &
   aap_tsq_k, r_bbp_tsq_k, ccp_tsq_k,                                          &
   ppp_tc_k(1, 0), ppp_tc_k(1, 1), ppp_tc_k(1, 2),                             &
   aap_qsq_k, r_bbp_qsq_k, ccp_qsq_k,                                          &
   ppp_qc_k(1, 0), ppp_qc_k(1, 1), ppp_qc_k(1, 2),                             &
   aap_cov_k, r_bbp_cov_k, ccp_cov_k,                                          &
   ppp_ct_k(1, 0), ppp_cq_k(1, 0),                                             &
   ppp_ct_k(1, 1), ppp_cq_k(1, 1),                                             &
   ppp_ct_k(1, 2), ppp_cq_k(1, 2))

r_qq_norm = 0.0
alp_num = 0.0
do k = 1, tke_levels
  ! set the initial values
  tsq_k(k) = 0.0
  qsq_k(k) = 0.0
  cov_k(k) = 0.0

  ! rvec is a residual vector
  rvec_tsq(k) = qq_tsq_k(k)
  rvec_qsq(k) = qq_qsq_k(k)
  rvec_cov(k) = qq_cov_k(k)

  r0vec_tsq(k) = rvec_tsq(k)
  r0vec_qsq(k) = rvec_qsq(k)
  r0vec_cov(k) = rvec_cov(k)

  pvec_tsq(k) = rvec_tsq(k)
  pvec_qsq(k) = rvec_qsq(k)
  pvec_cov(k) = rvec_cov(k)

  alp_num = alp_num + r0vec_tsq(k) * rvec_tsq(k)                               &
                    + r0vec_qsq(k) * rvec_qsq(k)                               &
                    + r0vec_cov(k) * rvec_cov(k)

  r_qq_norm = r_qq_norm +  qq_tsq_k(k) * qq_tsq_k(k)                           &
                        +  qq_qsq_k(k) * qq_qsq_k(k)                           &
                        +  qq_cov_k(k) * qq_cov_k(k)

end do

if (r_qq_norm == 0.0) then
  r_qq_norm = 0.0
  endflag = 2
  nitr = 0
else
  r_qq_norm = 1.0 / r_qq_norm
  endflag = -1
  nitr = max_itr
end if

outer_m_loop: do m = 1, nitr
  call mym_solve_simeq_ilud2(                                                  &
     0,                                                                        &
     pvec_tsq, pvec_qsq, pvec_cov,                                             &
     aap_tsq_k, r_bbp_tsq_k, ccp_tsq_k,                                        &
     ppp_tc_k(1, 0), ppp_tc_k(1, 1), ppp_tc_k(1, 2),                           &
     aap_qsq_k, r_bbp_qsq_k, ccp_qsq_k,                                        &
     ppp_qc_k(1, 0), ppp_qc_k(1, 1), ppp_qc_k(1, 2),                           &
     aap_cov_k, r_bbp_cov_k, ccp_cov_k,                                        &
     ppp_ct_k(1, 0), ppp_cq_k(1, 0),                                           &
     ppp_ct_k(1, 1), ppp_cq_k(1, 1),                                           &
     ppp_ct_k(1, 2), ppp_cq_k(1, 2),                                           &
     ppvec_tsq, ppvec_qsq, ppvec_cov)

  ! v = A pp
  call mym_simeq_matrix_prod(                                                  &
     aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                    &
     aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                    &
     aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                           &
     ppvec_tsq, ppvec_qsq, ppvec_cov,                                          &
     vvec_tsq, vvec_qsq, vvec_cov)

  alp_den = 0.0
  do k = 1, tke_levels
    alp_den = alp_den + r0vec_tsq(k) * vvec_tsq(k)                             &
       + r0vec_qsq(k) * vvec_qsq(k)                                            &
       + r0vec_cov(k) * vvec_cov(k)
  end do

  if (alp_den == 0.0) then
    endflag = 1
  else
    alp = alp_num / alp_den

    do k = 1, tke_levels
      svec_tsq(k) = rvec_tsq(k) - alp * vvec_tsq(k)
      svec_qsq(k) = rvec_qsq(k) - alp * vvec_qsq(k)
      svec_cov(k) = rvec_cov(k) - alp * vvec_cov(k)
    end do

    call mym_solve_simeq_ilud2(                                                &
       0,                                                                      &
       svec_tsq, svec_qsq, svec_cov,                                           &
       aap_tsq_k, r_bbp_tsq_k, ccp_tsq_k,                                      &
       ppp_tc_k(1, 0), ppp_tc_k(1, 1), ppp_tc_k(1, 2),                         &
       aap_qsq_k, r_bbp_qsq_k, ccp_qsq_k,                                      &
       ppp_qc_k(1, 0), ppp_qc_k(1, 1), ppp_qc_k(1, 2),                         &
       aap_cov_k, r_bbp_cov_k, ccp_cov_k,                                      &
       ppp_ct_k(1, 0), ppp_cq_k(1, 0),                                         &
       ppp_ct_k(1, 1), ppp_cq_k(1, 1),                                         &
       ppp_ct_k(1, 2), ppp_cq_k(1, 2),                                         &
       ssvec_tsq, ssvec_qsq, ssvec_cov)

    ! t = A ss
    call mym_simeq_matrix_prod(                                                &
       aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                  &
       aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                  &
       aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                         &
       ssvec_tsq, ssvec_qsq, ssvec_cov,                                        &
       tvec_tsq, tvec_qsq, tvec_cov)

    omg_num = 0.0
    omg_den = 0.0
    do k = 1, tke_levels
      omg_num = omg_num + tvec_tsq(k) * svec_tsq(k)                            &
                        + tvec_qsq(k) * svec_qsq(k)                            &
                        + tvec_cov(k) * svec_cov(k)
      omg_den = omg_den + tvec_tsq(k) * tvec_tsq(k)                            &
                        + tvec_qsq(k) * tvec_qsq(k)                            &
                        + tvec_cov(k) * tvec_cov(k)
    end do

    omg = omg_num / omg_den

    alp_den = alp_num

    alp_num = 0.0
    norm = 0.0
    max_val = 0.0
    do k = 1, tke_levels
      tsq_k(k) = tsq_k(k) + alp * ppvec_tsq(k) + omg * ssvec_tsq(k)
      qsq_k(k) = qsq_k(k) + alp * ppvec_qsq(k) + omg * ssvec_qsq(k)
      cov_k(k) = cov_k(k) + alp * ppvec_cov(k) + omg * ssvec_cov(k)
      rvec_tsq(k) = svec_tsq(k) - omg * tvec_tsq(k)
      rvec_qsq(k) = svec_qsq(k) - omg * tvec_qsq(k)
      rvec_cov(k) = svec_cov(k) - omg * tvec_cov(k)

      alp_num = alp_num + r0vec_tsq(k) * rvec_tsq(k)                           &
                        + r0vec_qsq(k) * rvec_qsq(k)                           &
                        + r0vec_cov(k) * rvec_cov(k)
      norm = norm + rvec_tsq(k) * rvec_tsq(k)                                  &
                  + rvec_qsq(k) * rvec_qsq(k)                                  &
                  + rvec_cov(k) * rvec_cov(k)

      max_val = max(max_val, abs(tsq_k(k)),                                    &
                             abs(qsq_k(k)),                                    &
                             abs(cov_k(k)))
    end do
    err = sqrt(norm * r_qq_norm)

    if (err >= eps .and. m < 30 .and. max_val < 1.0e10) then
      ! continue to the next step
    else if (max_val > 100.0) then
      ! Unexpectedly huge
      endflag = -2
    else if (err < eps) then
      ! Converged
      endflag = 0
    end if
  end if
  if (endflag /= -1) then
    exit outer_m_loop
  else
    bet = alp_num * alp / (alp_den * omg)
    do k = 1, tke_levels
      pvec_tsq(k) = rvec_tsq(k)                                                &
                           + bet * (pvec_tsq(k) - omg * vvec_tsq(k))
      pvec_qsq(k) = rvec_qsq(k)                                                &
                           + bet * (pvec_qsq(k) - omg * vvec_qsq(k))
      pvec_cov(k) = rvec_cov(k)                                                &
                           + bet * (pvec_cov(k) - omg * vvec_cov(k))
    end do
  end if
end do outer_m_loop    ! loop m = 1, max_itr

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_solve_simeq_bcgstab
end module mym_solve_simeq_bcgstab_mod
