!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Do gaussian elimination of hepta diagonal matrix.
module mym_simeq_solve_hepta_mod

  implicit none

  private
  public :: mym_simeq_solve_hepta

contains

  !> @brief Do gaussian elimination of hepta diagonal matrix.
  !> @param[in]     dd                  Diagonal element of matrix
  !> @param[in]     ee                  First off-diagonal element of matrix
  !> @param[in]     ff                  Second off-diagonal element of matrix
  !> @param[in]     gg                  Third off-diagonal element of matrix
  !> @param[in]     aa                  Third off-diagonal element of matrix
  !> @param[in]     bb                  Second off-diagonal element of matrix
  !> @param[in]     cc                  First off-diagonal element of matrix
  !> @param[in]     qq                  Right hand side term of equation
  subroutine mym_simeq_solve_hepta(dd, ee, ff, gg, aa, bb, cc, qq)

    use constants_mod, only: i_def, r_bl
    use bl_option_mod, only: one
    use mym_option_mod, only: tke_levels

    implicit none

    real(kind=r_bl), dimension(3*tke_levels), intent(in) :: dd, ee, ff, gg

    real(kind=r_bl), dimension(3*tke_levels), intent(inout) :: aa, bb, cc, qq

    integer(kind=i_def) :: k

    real(kind=r_bl) :: ddk, rddk, fac1, fac2, fac3

    k = 3*tke_levels
    rddk = one / dd(k)
    aa(k) = aa(k) * rddk
    bb(k) = bb(k) * rddk
    cc(k) = cc(k) * rddk
    qq(k) = qq(k) * rddk

    k = 3*tke_levels - 1
    fac1 = ee(k)
    ddk = dd(k) - fac1 * cc(k + 1)
    rddk = one / ddk
    aa(k) = aa(k) * rddk
    bb(k) = ( bb(k) - fac1 * aa(k + 1) ) * rddk
    cc(k) = ( cc(k) - fac1 * bb(k + 1) ) * rddk
    qq(k) = ( qq(k) - fac1 * qq(k + 1) ) * rddk

    k = 3*tke_levels - 2
    fac2 = ff(k)
    fac1 = ee(k) - fac2 * cc(k + 2)
    ddk = dd(k) - fac1 * cc(k + 1) - fac2 * bb(k + 2)
    rddk = one / ddk
    aa(k) = aa(k) * rddk
    bb(k) = ( bb(k) - fac1 * aa(k + 1) ) * rddk
    cc(k) = ( cc(k) - fac1 * bb(k + 1) - fac2 * aa(k + 2) ) * rddk
    qq(k) = ( qq(k) - fac1 * qq(k + 1) - fac2 * qq(k + 2) ) * rddk

    do k = 3*tke_levels - 3, 4, -1
      fac3 = gg(k)
      fac2 = ff(k) - fac3 * cc(k + 3)
      fac1 = ee(k) - fac2 * cc(k + 2) - fac3 * bb(k + 3)
      ddk = dd(k) - fac1 * cc(k + 1) - fac2 * bb(k + 2) - fac3 * aa(k + 3)
      rddk = one / ddk
      aa(k) = aa(k) * rddk
      bb(k) = ( bb(k) - fac1 * aa(k + 1) ) * rddk
      cc(k) = ( cc(k) - fac1 * bb(k + 1) - fac2 * aa(k + 2) ) * rddk
      qq(k) = ( qq(k) - fac1 * qq(k + 1) - fac2 * qq(k + 2)                   &
                                         - fac3 * qq(k + 3) ) * rddk
    end do

    k = 3
    fac3 = gg(k)
    fac2 = ff(k) - fac3 * cc(k + 3)
    fac1 = ee(k) - fac2 * cc(k + 2) - fac3 * bb(k + 3)
    ddk = dd(k) - fac1 * cc(k + 1) - fac2 * bb(k + 2) - fac3 * aa(k + 3)
    rddk = one / ddk
    bb(k) = ( bb(k) - fac1 * aa(k + 1) ) * rddk
    cc(k) = ( cc(k) - fac1 * bb(k + 1) - fac2 * aa(k + 2) ) * rddk
    qq(k) = ( qq(k) - fac1 * qq(k + 1) - fac2 * qq(k + 2)                     &
                                       - fac3 * qq(k + 3) ) * rddk

    k = 2
    fac3 = gg(k)
    fac2 = ff(k) - fac3 * cc(k + 3)
    fac1 = ee(k) - fac2 * cc(k + 2) - fac3 * bb(k + 3)
    ddk = dd(k) - fac1 * cc(k + 1) - fac2 * bb(k + 2) - fac3 * aa(k + 3)
    rddk = one / ddk
    cc(k) = ( cc(k) - fac1 * bb(k + 1) - fac2 * aa(k + 2) ) * rddk
    qq(k) = ( qq(k) - fac1 * qq(k + 1) - fac2 * qq(k + 2)                     &
                                       - fac3 * qq(k + 3) ) * rddk

    k = 1
    fac3 = gg(k)
    fac2 = ff(k) - fac3 * cc(k + 3)
    fac1 = ee(k) - fac2 * cc(k + 2) - fac3 * bb(k + 3)
    ddk = dd(k) - fac1 * cc(k + 1) - fac2 * bb(k + 2) - fac3 * aa(k + 3)
    rddk = one / ddk
    qq(k) = ( qq(k) - fac1 * qq(k + 1) - fac2 * qq(k + 2)                     &
                                       - fac3 * qq(k + 3) ) * rddk

    k = 2
    qq(k) = qq(k) - cc(k) * qq(k - 1)

    k = 3
    qq(k) = qq(k) - cc(k) * qq(k - 1) - bb(k) * qq(k - 2)

    do k = 4, 3*tke_levels
      qq(k) = qq(k) - cc(k) * qq(k - 1) - bb(k) * qq(k - 2) - aa(k) * qq(k - 3)
    end do

    return

  end subroutine mym_simeq_solve_hepta

end module mym_simeq_solve_hepta_mod
