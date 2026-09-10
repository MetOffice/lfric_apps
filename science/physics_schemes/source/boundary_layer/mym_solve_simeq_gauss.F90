!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Solve simultaneous equations with gaussian elimination.
module mym_solve_simeq_gauss_mod

  implicit none

  private
  public :: mym_solve_simeq_gauss

contains
  !> @brief Solve simultaneous equations with gaussian elimination.
  !> @param[in]     qq_tsq_k            Right hand side term of tsq(k) equation
  !> @param[in]     qq_qsq_k            Right hand side term of qsq(k) equation
  !> @param[in]     qq_cov_k            Right hand side term of cov(k) equation
  !> @param[in]     aa_tsq_k            Coefficient of tsq(k-1)
  !> @param[in]     bb_tsq_k            Coefficient of tsq(k)
  !> @param[in]     cc_tsq_k            Coefficient of tsq(k+1)
  !> @param[in]     pp_tc_k             Correlation with cov in tsq equation
  !> @param[in]     aa_qsq_k            Coefficient of qsq(k-1)
  !> @param[in]     bb_qsq_k            Coefficient of qsq(k)
  !> @param[in]     cc_qsq_k            Coefficient of qsq(k+1)
  !> @param[in]     pp_qc_k             Correlation with cov in qsq equation
  !> @param[in]     aa_cov_k            Coefficient of cov(k-1)
  !> @param[in]     bb_cov_k            Coefficient of cov(k)
  !> @param[in]     cc_cov_k            Coefficient of cov(k+1)
  !> @param[in]     pp_ct_k             Correlation with tsq in cov equation
  !> @param[in]     pp_cq_k             Correlation with qsq in cov equation
  !> @param[out]    tsq_k               Solved tsq
  !> @param[out]    qsq_k               Solved qsq
  !> @param[out]    cov_k               Solved cov
  subroutine mym_solve_simeq_gauss(                                            &
     qq_tsq_k, qq_qsq_k, qq_cov_k,                                             &
     aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                                    &
     aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                                    &
     aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                           &
     tsq_k, qsq_k, cov_k)

    use constants_mod, only: i_def, r_bl
    use bl_option_mod, only: zero
    use mym_option_mod, only: tke_levels
    use mym_simeq_solve_hepta_mod, only: mym_simeq_solve_hepta

    implicit none

    real(kind=r_bl), dimension(tke_levels), intent(in) :: qq_tsq_k,            &
                                                          qq_qsq_k,            &
                                                          qq_cov_k,            &
                                                          aa_tsq_k,            &
                                                          bb_tsq_k,            &
                                                          cc_tsq_k,            &
                                                          pp_tc_k,             &
                                                          aa_qsq_k,            &
                                                          bb_qsq_k,            &
                                                          cc_qsq_k,            &
                                                          pp_qc_k,             &
                                                          aa_cov_k,            &
                                                          bb_cov_k,            &
                                                          cc_cov_k,            &
                                                          pp_ct_k,             &
                                                          pp_cq_k

    real(kind=r_bl), dimension(tke_levels), intent(out) :: tsq_k,              &
                                                           qsq_k,              &
                                                           cov_k
    ! Local variables
    integer(kind=i_def) :: k, kt, kq, kc

    real(kind=r_bl), dimension(3*tke_levels) :: aa, bb, cc, dd, ee, ff, gg, qq
             ! elements of heptadiag matrix

    k = 1
    kt = 1
    kq = 2
    kc = 3

    bb(kc) = pp_ct_k(k)

    cc(kq) = zero
    cc(kc) = pp_cq_k(k)

    dd(kt) = bb_tsq_k(k)
    dd(kq) = bb_qsq_k(k)
    dd(kc) = bb_cov_k(k)

    ee(kt) = zero
    ee(kq) = pp_qc_k(k)
    ee(kc) = zero

    ff(kt) = pp_tc_k(k)
    ff(kq) = zero
    ff(kc) = zero

    gg(kt) = cc_tsq_k(k)
    gg(kq) = cc_qsq_k(k)
    gg(kc) = cc_cov_k(k)

    qq(kt) = qq_tsq_k(k)
    qq(kq) = qq_qsq_k(k)
    qq(kc) = qq_cov_k(k)

    do k = 2, tke_levels - 1
      kt = 3*k - 2
      kq = kt + 1
      kc = kt + 2

      aa(kt) = aa_tsq_k(k)
      aa(kq) = aa_qsq_k(k)
      aa(kc) = aa_cov_k(k)

      bb(kt) = zero
      bb(kq) = zero
      bb(kc) = pp_ct_k(k)

      cc(kt) = zero
      cc(kq) = zero
      cc(kc) = pp_cq_k(k)

      dd(kt) = bb_tsq_k(k)
      dd(kq) = bb_qsq_k(k)
      dd(kc) = bb_cov_k(k)

      ee(kt) = zero
      ee(kq) = pp_qc_k(k)
      ee(kc) = zero

      ff(kt) = pp_tc_k(k)
      ff(kq) = zero
      ff(kc) = zero

      gg(kt) = cc_tsq_k(k)
      gg(kq) = cc_qsq_k(k)
      gg(kc) = cc_cov_k(k)

      qq(kt) = qq_tsq_k(k)
      qq(kq) = qq_qsq_k(k)
      qq(kc) = qq_cov_k(k)
    end do

    k = tke_levels
    kt = 3*k - 2
    kq = kt + 1
    kc = kt + 2

    aa(kt) = aa_tsq_k(k)
    aa(kq) = aa_qsq_k(k)
    aa(kc) = aa_cov_k(k)

    bb(kt) = zero
    bb(kq) = zero
    bb(kc) = pp_ct_k(k)

    cc(kt) = zero
    cc(kq) = zero
    cc(kc) = pp_cq_k(k)

    dd(kt) = bb_tsq_k(k)
    dd(kq) = bb_qsq_k(k)
    dd(kc) = bb_cov_k(k)

    ee(kt) = zero
    ee(kq) = pp_qc_k(k)

    ff(kt) = pp_tc_k(k)

    qq(kt) = qq_tsq_k(k)
    qq(kq) = qq_qsq_k(k)
    qq(kc) = qq_cov_k(k)

    call mym_simeq_solve_hepta(dd, ee, ff, gg, aa, bb, cc, qq)

    do k = 1, tke_levels
      kt = 3*k - 2
      kq = kt + 1
      kc = kt + 2

      tsq_k(k) = qq(kt)
      qsq_k(k) = qq(kq)
      cov_k(k) = qq(kc)
    end do

    return

  end subroutine mym_solve_simeq_gauss

end module mym_solve_simeq_gauss_mod
