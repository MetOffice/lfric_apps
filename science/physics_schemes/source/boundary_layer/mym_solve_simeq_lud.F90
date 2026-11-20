! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To solve simultaneous equations by LU decomposition

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_solve_simeq_lud_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_SOLVE_SIMEQ_LUD_MOD'
contains

subroutine mym_solve_simeq_lud(                                                &
      qq_tsq_k, qq_qsq_k, qq_cov_k,                                            &
      aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                   &
      aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                   &
      aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                          &
      tsq_k, qsq_k, cov_k)

use mym_option_mod, only: tke_levels
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none

! intent in variables
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

integer ::                                                                     &
   k, l, m, n,                                                                 &
           ! loop indexes
   kpiv
           ! index of a pivot

real(kind=r_bl) ::                                                      &
   wk
          ! work variables

real(kind=r_bl) ::                                                      &
   amat(3 * tke_levels, 3 * tke_levels),                                       &
           ! coefficient matrix
   bvec(3 * tke_levels)
           ! vector in the right hand side

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_SOLVE_SIMEQ_LUD'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

amat(:, :) = 0.0

do k = 1, tke_levels
  amat(k, k) = bb_tsq_k(k)
  amat(tke_levels + k, tke_levels + k) = bb_qsq_k(k)
  amat(2 * tke_levels + k, 2 * tke_levels + k)                                 &
                                        = bb_cov_k(k)
  bvec(k)    = qq_tsq_k(k)
  bvec(tke_levels + k) = qq_qsq_k(k)
  bvec(2 * tke_levels + k) = qq_cov_k(k)
end do

do k = 2, tke_levels
  amat(k, k-1) = aa_tsq_k(k)
  amat(tke_levels + k, tke_levels + k - 1) = aa_qsq_k(k)
  amat(2 * tke_levels + k, 2 * tke_levels + k - 1)                             &
                                        = aa_cov_k(k)
end do

do k = 1, tke_levels - 1
  amat(k, k+1) = cc_tsq_k(k)
  amat(tke_levels + k, tke_levels + k + 1) = cc_qsq_k(k)
  amat(2 * tke_levels + k, 2 * tke_levels + k + 1)                             &
                                        = cc_cov_k(k)
end do

do k = 1, tke_levels
  amat(k, 2 * tke_levels + k) = pp_tc_k(k)
  amat(tke_levels + k, 2 * tke_levels + k) = pp_qc_k(k)
  amat(2 * tke_levels + k, k) = pp_ct_k(k)
  amat(2 * tke_levels + k, tke_levels + k) = pp_cq_k(k)
end do

n = 3 * tke_levels
! main part
do k = 1, n
  kpiv = k
  wk  = abs(amat(k, k))
  do l = k + 1, n
    if (abs(amat(l, k)) > wk) then
      kpiv = l
      wk  = abs(amat(l, k))
    end if
  end do

  if (kpiv /= k) then
    do m = 1, n
      wk       = amat(k, m)
      amat(k, m)    = amat(kpiv, m)
      amat(kpiv, m) = wk
    end do
    wk   = bvec(k)
    bvec(k) = bvec(kpiv)
    bvec(kpiv) = wk
  end if

  amat(k, k) = 1.0 / amat(k, k)

  do l = k + 1, n
    amat(l, k) = amat(l, k) * amat(k, k)
  end do

  do m = k + 1, n
    do l = k+1, n
      amat(l, m) = amat(l, m) - amat(k, m) * amat(l, k)
    end do
  end do
end do  ! loop k = 1, n

do m = 1, n - 1
  do l = m + 1, n
    bvec(l) = bvec(l) - bvec(m) * amat(l, m)
  end do
end do

do m = n, 1, -1
  bvec(m) = bvec(m) * amat(m, m)
  do l = 1, m - 1
    bvec(l) = bvec(l) - amat(l, m) * bvec(m)
  end do
end do

do k = 1, tke_levels
  tsq_k(k) = bvec(k)
  qsq_k(k) = bvec(tke_levels + k)
  cov_k(k) = bvec(2 * tke_levels + k)
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_solve_simeq_lud
end module mym_solve_simeq_lud_mod
