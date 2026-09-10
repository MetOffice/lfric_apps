! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To perform the incomplete LU decomposition with fill-in
!           level 2

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_simeq_ilud2_decmp_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_SIMEQ_ILUD2_DECMP_MOD'
contains

subroutine mym_simeq_ilud2_decmp(                                              &
      aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                   &
      aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                   &
      aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                          &
      aap_tsq_k, r_bbp_tsq_k, ccp_tsq_k,                                       &
      ppp_tc_k, pp1_tc_k,  pp2_tc_k,                                           &
      aap_qsq_k, r_bbp_qsq_k, ccp_qsq_k,                                       &
      ppp_qc_k, pp1_qc_k,  pp2_qc_k,                                           &
      aap_cov_k, r_bbp_cov_k, ccp_cov_k,                                       &
      ppp_ct_k, ppp_cq_k, pp1_ct_k, pp1_cq_k, pp2_ct_k, pp2_cq_k)

use mym_option_mod, only: tke_levels
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none

! intent in variables
real(kind=r_bl), intent(in) ::                                          &
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
           ! matrix elements of the ILU decomposed matrix

integer :: k
           ! loop indexes

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_SIMEQ_ILUD2_DECMP'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

aap_tsq_k(1) = aa_tsq_k(1)
r_bbp_tsq_k(1) = 1.0 / bb_tsq_k(1)
ccp_tsq_k(1) = cc_tsq_k(1)

aap_qsq_k(1) = aa_qsq_k(1)
r_bbp_qsq_k(1) = 1.0 / bb_qsq_k(1)
ccp_qsq_k(1) = cc_qsq_k(1)

ppp_tc_k(1) = pp_tc_k(1)
ppp_qc_k(1) = pp_qc_k(1)

pp1_tc_k(1) = 0.0
pp1_qc_k(1) = 0.0
pp2_tc_k(1) = 0.0
pp2_qc_k(1) = 0.0

do k = 2, tke_levels
  aap_tsq_k(k) = aa_tsq_k(k)
  r_bbp_tsq_k(k) = 1.0 / (bb_tsq_k(k)                                          &
            - aap_tsq_k(k) * ccp_tsq_k(k - 1) * r_bbp_tsq_k(k - 1))
  ccp_tsq_k(k) = cc_tsq_k(k)

  aap_qsq_k(k) = aa_qsq_k(k)
  r_bbp_qsq_k(k) = 1.0 / (bb_qsq_k(k)                                          &
            - aap_qsq_k(k) * ccp_qsq_k(k - 1) * r_bbp_qsq_k(k - 1))
  ccp_qsq_k(k) = cc_qsq_k(k)

  ppp_tc_k(k) = pp_tc_k(k)
  ppp_qc_k(k) = pp_qc_k(k)

  pp1_tc_k(k) = - aap_tsq_k(k) * r_bbp_tsq_k(k - 1) * ppp_tc_k(k - 1)
  pp1_qc_k(k) = - aap_qsq_k(k) * r_bbp_qsq_k(k - 1) * ppp_qc_k(k - 1)

  pp2_tc_k(k) = - aap_tsq_k(k) * r_bbp_tsq_k(k - 1) * pp1_tc_k(k - 1)
  pp2_qc_k(k) = - aap_qsq_k(k) * r_bbp_qsq_k(k - 1) * pp1_qc_k(k - 1)

end do

k = 1

ppp_ct_k(k) = pp_ct_k(k)
ppp_cq_k(k) = pp_cq_k(k)

pp1_ct_k(k) = -ppp_ct_k(k) * r_bbp_tsq_k(k) * ccp_tsq_k(k)
pp1_cq_k(k) = -ppp_cq_k(k) * r_bbp_qsq_k(k) * ccp_qsq_k(k)

pp2_ct_k(k) = - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ccp_tsq_k(k + 1)
pp2_cq_k(k) = - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ccp_qsq_k(k + 1)

aap_cov_k(k) = 0.0

r_bbp_cov_k(k) = 1.0 / (                                                       &
                bb_cov_k(k)                                                    &
                - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)                   &
                - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp1_tc_k(k + 1)           &
                - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp2_tc_k(k + 2)           &
                - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)                   &
                - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp1_qc_k(k + 1)           &
                - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp2_qc_k(k + 2))

ccp_cov_k(k) = cc_cov_k(k)                                                     &
               - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ppp_tc_k(k + 1)            &
               - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp1_tc_k(k + 2)            &
               - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ppp_qc_k(k + 1)            &
               - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp1_qc_k(k + 2)

do k = 2, tke_levels - 2
  ppp_ct_k(k) = pp_ct_k(k)
  ppp_cq_k(k) = pp_cq_k(k)

  pp1_ct_k(k) = -ppp_ct_k(k) * r_bbp_tsq_k(k) * ccp_tsq_k(k)
  pp1_cq_k(k) = -ppp_cq_k(k) * r_bbp_qsq_k(k) * ccp_qsq_k(k)

  pp2_ct_k(k) = - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ccp_tsq_k(k + 1)
  pp2_cq_k(k) = - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ccp_qsq_k(k + 1)

  aap_cov_k(k) = aa_cov_k(k)                                                   &
                 - ppp_ct_k(k) * r_bbp_tsq_k(k) * pp1_tc_k(k)                  &
                 - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp2_tc_k(k + 1)          &
                 - ppp_cq_k(k) * r_bbp_qsq_k(k) * pp1_qc_k(k)                  &
                 - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp2_qc_k(k + 1)

  r_bbp_cov_k(k) = 1.0 / (                                                     &
               bb_cov_k(k)                                                     &
               - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)                    &
               - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp1_tc_k(k + 1)            &
               - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp2_tc_k(k + 2)            &
               - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)                    &
               - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp1_qc_k(k + 1)            &
               - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp2_qc_k(k + 2)            &
               - aap_cov_k(k) * r_bbp_cov_k(k - 1) * ccp_cov_k(k - 1))

  ccp_cov_k(k) = cc_cov_k(k)                                                   &
                 - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ppp_tc_k(k + 1)          &
                 - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp1_tc_k(k + 2)          &
                 - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ppp_qc_k(k + 1)          &
                 - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp1_qc_k(k + 2)
end do

k = tke_levels - 1

ppp_ct_k(k) = pp_ct_k(k)
ppp_cq_k(k) = pp_cq_k(k)

pp1_ct_k(k) = -ppp_ct_k(k) * r_bbp_tsq_k(k) * ccp_tsq_k(k)
pp1_cq_k(k) = -ppp_cq_k(k) * r_bbp_qsq_k(k) * ccp_qsq_k(k)
pp2_ct_k(k) = 0.0
pp2_cq_k(k) = 0.0

aap_cov_k(k) = aa_cov_k(k)                                                     &
               - ppp_ct_k(k) * r_bbp_tsq_k(k) * pp1_tc_k(k)                    &
               - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp2_tc_k(k + 1)            &
               - ppp_cq_k(k) * r_bbp_qsq_k(k) * pp1_qc_k(k)                    &
               - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp2_qc_k(k + 1)

r_bbp_cov_k(k) = 1.0 / (                                                       &
               bb_cov_k(k)                                                     &
               - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)                    &
               - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp1_tc_k(k + 1)            &
               - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)                    &
               - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp1_qc_k(k + 1)            &
               - aap_cov_k(k) * r_bbp_cov_k(k - 1) * ccp_cov_k(k - 1))

ccp_cov_k(k) = cc_cov_k(k)                                                     &
               - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ppp_tc_k(k + 1)            &
               - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ppp_qc_k(k + 1)

k = tke_levels

ppp_ct_k(k) = pp_ct_k(k)
ppp_cq_k(k) = pp_cq_k(k)

pp1_ct_k(k) = 0.0
pp1_cq_k(k) = 0.0
pp2_ct_k(k) = 0.0
pp2_cq_k(k) = 0.0

aap_cov_k(k) = aa_cov_k(k)                                                     &
               - ppp_ct_k(k) * r_bbp_tsq_k(k) * pp1_tc_k(k)                    &
               - ppp_cq_k(k) * r_bbp_qsq_k(k) * pp1_qc_k(k)
r_bbp_cov_k(k) = 1.0 / (                                                       &
               bb_cov_k(k)                                                     &
               - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)                    &
               - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)                    &
               - aap_cov_k(k) * r_bbp_cov_k(k - 1) * ccp_cov_k(k - 1))
ccp_cov_k(k) = 0.0

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_simeq_ilud2_decmp
end module mym_simeq_ilud2_decmp_mod
