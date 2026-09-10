! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate products of a matrix and a vector
!           in solving the simultaneous equations.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_simeq_matrix_prod_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_SIMEQ_MATRIX_PROD_MOD'
contains

subroutine mym_simeq_matrix_prod(                                              &
      aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                   &
      aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                   &
      aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                          &
      x_tsq_k, x_qsq_k, x_cov_k,                                               &
      y_tsq_k, y_qsq_k, y_cov_k)
use mym_option_mod, only: tke_levels
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none

real(kind=r_bl), intent(in) ::                                          &
   ! matrix elements (for meanings of each, see the document)
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
   ! vector elements
   x_tsq_k(tke_levels),                                                        &
   x_qsq_k(tke_levels),                                                        &
   x_cov_k(tke_levels)

real(kind=r_bl), intent(out) ::                                         &
   ! vector elements of products (answers)
   y_tsq_k(tke_levels),                                                        &
   y_qsq_k(tke_levels),                                                        &
   y_cov_k(tke_levels)

integer :: k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_SIMEQ_MATRIX_PROD'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! y = A * x
k = 1
y_tsq_k(k) =      bb_tsq_k(k) * x_tsq_k(k)                                     &
                + cc_tsq_k(k) * x_tsq_k(k + 1)                                 &
                + pp_tc_k(k) * x_cov_k(k)

y_qsq_k(k) =      bb_qsq_k(k) * x_qsq_k(k)                                     &
                + cc_qsq_k(k) * x_qsq_k(k + 1)                                 &
                + pp_qc_k(k) * x_cov_k(k)

y_cov_k(k) =     bb_cov_k(k) * x_cov_k(k)                                      &
                + cc_cov_k(k) * x_cov_k(k + 1)                                 &
                + pp_ct_k(k) * x_tsq_k(k)                                      &
                + pp_cq_k(k) * x_qsq_k(k)

do k = 2, tke_levels - 1
  y_tsq_k(k) = aa_tsq_k(k) * x_tsq_k(k - 1)                                    &
                + bb_tsq_k(k) * x_tsq_k(k)                                     &
                + cc_tsq_k(k) * x_tsq_k(k + 1)                                 &
                + pp_tc_k(k) * x_cov_k(k)

  y_qsq_k(k) = aa_qsq_k(k) * x_qsq_k(k - 1)                                    &
                + bb_qsq_k(k) * x_qsq_k(k)                                     &
                + cc_qsq_k(k) * x_qsq_k(k + 1)                                 &
                + pp_qc_k(k) * x_cov_k(k)

  y_cov_k(k) = aa_cov_k(k) * x_cov_k(k - 1)                                    &
                + bb_cov_k(k) * x_cov_k(k)                                     &
                + cc_cov_k(k) * x_cov_k(k + 1)                                 &
                + pp_ct_k(k) * x_tsq_k(k)                                      &
                + pp_cq_k(k) * x_qsq_k(k)

end do

k = tke_levels
y_tsq_k(k) = aa_tsq_k(k) * x_tsq_k(k - 1)                                      &
                + bb_tsq_k(k) * x_tsq_k(k)                                     &
                + pp_tc_k(k) * x_cov_k(k)

y_qsq_k(k) = aa_qsq_k(k) * x_qsq_k(k - 1)                                      &
                + bb_qsq_k(k) * x_qsq_k(k)                                     &
                + pp_qc_k(k) * x_cov_k(k)

y_cov_k(k) = aa_cov_k(k) * x_cov_k(k - 1)                                      &
                + bb_cov_k(k) * x_cov_k(k)                                     &
                + pp_ct_k(k) * x_tsq_k(k)                                      &
                + pp_cq_k(k) * x_qsq_k(k)


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_simeq_matrix_prod
end module mym_simeq_matrix_prod_mod
