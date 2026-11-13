! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To solve simultaneous equations of which the coefficient
!           matrix is obtained by imcompelete LU decomposition
!           with fill-in level 2 (ILU(2)) for the original coefficient
!           matrix.
!           ILU(2) decomposition is assumed to have been already done
!           in mym_simeq_ilud2_dcmp.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_solve_simeq_ilud2_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_SOLVE_SIMEQ_ILUD2_MOD'
contains

subroutine mym_solve_simeq_ilud2(                                              &
      imode,                                                                   &
      qq_tsq_k, qq_qsq_k, qq_cov_k,                                            &
      aap_tsq_k, r_bbp_tsq_k, ccp_tsq_k,                                       &
      ppp_tc_k, pp1_tc_k,  pp2_tc_k,                                           &
      aap_qsq_k, r_bbp_qsq_k, ccp_qsq_k,                                       &
      ppp_qc_k, pp1_qc_k,  pp2_qc_k,                                           &
      aap_cov_k, r_bbp_cov_k, ccp_cov_k,                                       &
      ppp_ct_k, ppp_cq_k, pp1_ct_k, pp1_cq_k, pp2_ct_k, pp2_cq_k,              &
      tsq_k, qsq_k, cov_k)

use mym_option_mod, only: tke_levels
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none

! intent in variables
integer, intent(in) :: imode
                       ! mode switch for the Matrix
                       ! 0: normal, 1: transposed

real(kind=r_bl), intent(in) ::                                          &
   qq_tsq_k(tke_levels),                                                       &
   qq_qsq_k(tke_levels),                                                       &
   qq_cov_k(tke_levels),                                                       &
   aap_tsq_k(tke_levels),                                                      &
   r_bbp_tsq_k(tke_levels),                                                    &
   ccp_tsq_k(tke_levels),                                                      &
   ppp_tc_k(tke_levels),                                                       &
   pp1_tc_k(tke_levels),                                                       &
   pp2_tc_k(tke_levels),                                                       &
   aap_qsq_k(tke_levels),                                                      &
   r_bbp_qsq_k(tke_levels),                                                    &
   ccp_qsq_k(tke_levels),                                                      &
   ppp_qc_k(tke_levels),                                                       &
   pp1_qc_k(tke_levels),                                                       &
   pp2_qc_k(tke_levels),                                                       &
   aap_cov_k(tke_levels),                                                      &
   r_bbp_cov_k(tke_levels),                                                    &
   ccp_cov_k(tke_levels),                                                      &
   ppp_ct_k(tke_levels),                                                       &
   ppp_cq_k(tke_levels),                                                       &
   pp1_ct_k(tke_levels),                                                       &
   pp1_cq_k(tke_levels),                                                       &
   pp2_ct_k(tke_levels),                                                       &
   pp2_cq_k(tke_levels)
           ! matrix elements of ILU decomposed matrix
           ! See the document for details

real(kind=r_bl), intent(out) ::                                         &
   tsq_k(tke_levels),                                                          &
   qsq_k(tke_levels),                                                          &
   cov_k(tke_levels)
           ! solution vectors

integer :: k
           ! loop indexes

integer, parameter ::                                                          &
   normal = 0,                                                                 &
   transposed = 1
           ! symbols for the mode

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_SOLVE_SIMEQ_ILUD2'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tsq_k(1) = qq_tsq_k(1) * r_bbp_tsq_k(1)
qsq_k(1) = qq_qsq_k(1) * r_bbp_qsq_k(1)

if (imode == normal) then
  do k = 2, tke_levels
    tsq_k(k) = (qq_tsq_k(k) - aap_tsq_k(k) * tsq_k(k - 1))                     &
                                               * r_bbp_tsq_k(k)
    qsq_k(k) = (qq_qsq_k(k) - aap_qsq_k(k) * qsq_k(k - 1))                     &
                                               * r_bbp_qsq_k(k)
  end do

  k = 1
  cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                             &
                          - pp1_ct_k(k) * tsq_k(k + 1)                         &
                          - pp2_ct_k(k) * tsq_k(k + 2)                         &
                          - ppp_cq_k(k) * qsq_k(k)                             &
                          - pp1_cq_k(k) * qsq_k(k + 1)                         &
                          - pp2_cq_k(k) * qsq_k(k + 2))                        &
             * r_bbp_cov_k(k)

  do k = 2, tke_levels - 2
    cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                           &
                            - pp1_ct_k(k) * tsq_k(k + 1)                       &
                            - pp2_ct_k(k) * tsq_k(k + 2)                       &
                            - ppp_cq_k(k) * qsq_k(k)                           &
                            - pp1_cq_k(k) * qsq_k(k + 1)                       &
                            - pp2_cq_k(k) * qsq_k(k + 2)                       &
                            - aap_cov_k(k) * cov_k(k - 1))                     &
             * r_bbp_cov_k(k)
  end do

  k = tke_levels - 1
  cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                             &
                          - pp1_ct_k(k) * tsq_k(k + 1)                         &
                          - ppp_cq_k(k) * qsq_k(k)                             &
                          - pp1_cq_k(k) * qsq_k(k + 1)                         &
                          - aap_cov_k(k) * cov_k(k - 1))                       &
             * r_bbp_cov_k(k)


  k = tke_levels
  cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                             &
                          - ppp_cq_k(k) * qsq_k(k)                             &
                          - aap_cov_k(k) * cov_k(k - 1))                       &
             * r_bbp_cov_k(k)


  do k = tke_levels - 1, 1, -1
    cov_k(k) = cov_k(k)                                                        &
                      - ccp_cov_k(k) * cov_k(k + 1) * r_bbp_cov_k(k)
  end do

  k = tke_levels
  qsq_k(k) = qsq_k(k) - (ppp_qc_k(k) * cov_k(k)                                &
                        + pp1_qc_k(k) * cov_k(k - 1)                           &
                        + pp2_qc_k(k) * cov_k(k - 2))                          &
                       * r_bbp_qsq_k(k)
  tsq_k(k) = tsq_k(k) - (ppp_tc_k(k) * cov_k(k)                                &
                        + pp1_tc_k(k) * cov_k(k - 1)                           &
                        + pp2_tc_k(k) * cov_k(k - 2))                          &
                       * r_bbp_tsq_k(k)

  do k = tke_levels - 1, 3, -1
    qsq_k(k) = qsq_k(k)                                                        &
          - (ppp_qc_k(k) * cov_k(k) + ccp_qsq_k(k) * qsq_k(k + 1)              &
             + pp1_qc_k(k) * cov_k(k - 1)                                      &
             + pp2_qc_k(k) * cov_k(k - 2))                                     &
          * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k)                                                        &
          - (ppp_tc_k(k) * cov_k(k) + ccp_tsq_k(k) * tsq_k(k + 1)              &
             + pp1_tc_k(k) * cov_k(k - 1)                                      &
             + pp2_tc_k(k) * cov_k(k - 2))                                     &
          * r_bbp_tsq_k(k)
  end do

  k = 2
  qsq_k(k) = qsq_k(k)                                                          &
        - (ppp_qc_k(k) * cov_k(k) + ccp_qsq_k(k) * qsq_k(k + 1)                &
           + pp1_qc_k(k) * cov_k(k - 1))                                       &
        * r_bbp_qsq_k(k)
  tsq_k(k) = tsq_k(k)                                                          &
        - (ppp_tc_k(k) * cov_k(k) + ccp_tsq_k(k) * tsq_k(k + 1)                &
           + pp1_tc_k(k) * cov_k(k - 1))                                       &
        * r_bbp_tsq_k(k)


  k = 1
  qsq_k(k) = qsq_k(k)                                                          &
        - (ppp_qc_k(k) * cov_k(k) + ccp_qsq_k(k) * qsq_k(k + 1))               &
        * r_bbp_qsq_k(k)
  tsq_k(k) = tsq_k(k)                                                          &
        - (ppp_tc_k(k) * cov_k(k) + ccp_tsq_k(k) * tsq_k(k + 1))               &
        * r_bbp_tsq_k(k)

else if (imode == transposed) then
  do k = 2, tke_levels
    tsq_k(k) = (qq_tsq_k(k)                                                    &
                 - ccp_tsq_k(k - 1) * tsq_k(k - 1)) * r_bbp_tsq_k(k)
    qsq_k(k) = (qq_qsq_k(k)                                                    &
                 - ccp_qsq_k(k - 1) * qsq_k(k - 1)) * r_bbp_qsq_k(k)
  end do

  k = 1
  cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                             &
                          - pp1_tc_k(k + 1) * tsq_k(k + 1)                     &
                          - pp2_tc_k(k + 2) * tsq_k(k + 2)                     &
                          - ppp_qc_k(k) * qsq_k(k)                             &
                          - pp1_qc_k(k + 1) * qsq_k(k + 1)                     &
                          - pp2_qc_k(k + 2) * qsq_k(k + 2))                    &
            * r_bbp_cov_k(k)

  do k = 2, tke_levels - 2
    cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                           &
                            - pp1_tc_k(k + 1) * tsq_k(k + 1)                   &
                            - pp2_tc_k(k + 2) * tsq_k(k + 2)                   &
                            - ppp_qc_k(k) * qsq_k(k)                           &
                            - pp1_qc_k(k + 1) * qsq_k(k + 1)                   &
                            - pp2_qc_k(k + 2) * qsq_k(k + 2)                   &
                            - ccp_cov_k(k - 1) * cov_k(k - 1))                 &
              * r_bbp_cov_k(k)
  end do

  k = tke_levels - 1
  cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                             &
                          - pp1_tc_k(k + 1) * tsq_k(k + 1)                     &
                          - ppp_qc_k(k) * qsq_k(k)                             &
                          - pp1_qc_k(k + 1) * qsq_k(k + 1)                     &
                          - ccp_cov_k(k - 1) * cov_k(k - 1))                   &
            * r_bbp_cov_k(k)


  k = tke_levels
  cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                             &
                        - ppp_qc_k(k) * qsq_k(k)                               &
                        - ccp_cov_k(k - 1) * cov_k(k - 1))                     &
            * r_bbp_cov_k(k)

  do k = tke_levels - 1, 1, -1
    cov_k(k) = cov_k(k)                                                        &
                 - aap_cov_k(k + 1) * cov_k(k + 1) * r_bbp_cov_k(k)
  end do

  k = tke_levels
  qsq_k(k) = qsq_k(k) - (ppp_cq_k(k) * cov_k(k)                                &
                        + pp1_cq_k(k - 1) * cov_k(k - 1)                       &
                        + pp2_cq_k(k - 2) * cov_k(k - 2))                      &
                       * r_bbp_qsq_k(k)
  tsq_k(k) = tsq_k(k) - (ppp_ct_k(k) * cov_k(k)                                &
                        + pp1_ct_k(k - 1) * cov_k(k - 1)                       &
                        + pp2_ct_k(k - 2) * cov_k(k - 2))                      &
                       * r_bbp_tsq_k(k)

  do k = tke_levels - 1, 3, -1
    qsq_k(k) = qsq_k(k)                                                        &
          - (ppp_cq_k(k) * cov_k(k) + aap_qsq_k(k + 1) * qsq_k(k + 1)          &
             + pp1_cq_k(k - 1) * cov_k(k - 1)                                  &
             + pp2_cq_k(k - 2) * cov_k(k - 2))                                 &
          * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k)                                                        &
          - (ppp_ct_k(k) * cov_k(k) + aap_tsq_k(k + 1) * tsq_k(k + 1)          &
             + pp1_ct_k(k - 1) * cov_k(k - 1)                                  &
             + pp2_ct_k(k - 2) * cov_k(k - 2))                                 &
          * r_bbp_tsq_k(k)
  end do

  k = 2
  qsq_k(k) = qsq_k(k)                                                          &
        - (ppp_cq_k(k) * cov_k(k) + aap_qsq_k(k + 1) * qsq_k(k + 1)            &
           + pp1_cq_k(k - 1) * cov_k(k - 1))                                   &
        * r_bbp_qsq_k(k)
  tsq_k(k) = tsq_k(k)                                                          &
        - (ppp_ct_k(k) * cov_k(k) + aap_tsq_k(k + 1) * tsq_k(k + 1)            &
           + pp1_ct_k(k - 1) * cov_k(k - 1))                                   &
        * r_bbp_tsq_k(k)

  k = 1
  qsq_k(k) = qsq_k(k)                                                          &
        - (ppp_cq_k(k) * cov_k(k) + aap_qsq_k(k + 1) * qsq_k(k + 1))           &
        * r_bbp_qsq_k(k)
  tsq_k(k) = tsq_k(k)                                                          &
        - (ppp_ct_k(k) * cov_k(k) + aap_tsq_k(k + 1) * tsq_k(k + 1))           &
        * r_bbp_tsq_k(k)

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_solve_simeq_ilud2
end module mym_solve_simeq_ilud2_mod
